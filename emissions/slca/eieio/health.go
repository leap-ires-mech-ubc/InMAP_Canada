/*
Copyright © 2017 the InMAP authors.
This file is part of InMAP.

InMAP is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

InMAP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with InMAP.  If not, see <http://www.gnu.org/licenses/>.*/

package eieio

import (
	"context"
	"fmt"

	"github.com/ctessum/sparse"

	"github.com/spatialmodel/inmap/emissions/slca/eieio/eieiorpc"
	"github.com/spatialmodel/inmap/epi"
	"github.com/spatialmodel/inmap/internal/hash"

	"gonum.org/v1/gonum/mat"
)

type healthRequest struct {
	Demand     *mat.VecDense
	Industries *Mask
	Pol        Pollutant
	Pop        string
	Year       Year
	Loc        Location
	HR         epi.HRer
	AQM        string
}

// Health returns spatially-explicit pollutant air quality-related health impacts caused by the
// specified economic demand.  Emitters
// specify the emitters health impacts should be calculated for.
// If emitters == nil, combined health impacts for all emitters are calculated.
// Population must be one of the population types defined in the configuration file.
func (e *SpatialEIO) Health(ctx context.Context, request *eieiorpc.HealthInput) (*eieiorpc.Vector, error) {
	e.loadHealthOnce.Do(func() {
		var c string
		if e.EIEIOCache != "" {
			c = e.EIEIOCache + "/individual"
		}
		e.healthCache = loadCacheOnce(func(ctx context.Context, request interface{}) (interface{}, error) {
			r := request.(*healthRequest)
			return e.health(ctx, r.Demand, r.Industries, r.AQM, r.Pol, r.Pop, r.Year, r.Loc, r.HR) // Actually calculate the health impacts.
		}, 1, e.MemCacheSize, c, vectorMarshal, vectorUnmarshal)
	})
	hr, ok := e.hr[request.HR]
	if !ok {
		return nil, fmt.Errorf("eieio: hazard ratio function `%s` is not registered", request.HR)
	}
	req := &healthRequest{
		Demand:     array2vec(request.Demand.Data),
		Industries: rpc2mask(request.EmitterMask),
		Pol:        Pollutant(request.Pollutant),
		Pop:        request.Population,
		Year:       Year(request.Year),
		Loc:        Location(request.Location),
		HR:         hr,
		AQM:        request.AQM,
	}
	rr := e.healthCache.NewRequest(ctx, req, "health_"+hash.Hash(req))
	resultI, err := rr.Result()
	if err != nil {
		return nil, err
	}
	return vec2rpc(resultI.(*mat.VecDense)), nil
}

// health returns spatially-explicit pollutant air quality-related health impacts caused by the
// specified economic demand. Emitters
// specify the emitters health impacts should be calculated for.
// If emitters == nil, combined health impacts for all emitters are calculated.
// Population must be one of the population types defined in the configuration file.
func (e *SpatialEIO) health(ctx context.Context, demand *mat.VecDense, industries *Mask, aqm string, pol Pollutant, pop string, year Year, loc Location, HR epi.HRer) (*mat.VecDense, error) {
	hf, err := e.healthFactors(ctx, aqm, pol, pop, year, HR)
	if err != nil {
		return nil, err
	}

	activity, err := e.economicImpactsSCC(demand, year, loc)
	if err != nil {
		return nil, err
	}

	if industries != nil {
		// Set activity in industries we're not interested in to zero.
		industries.Mask(activity)
	}

	r, _ := hf.Dims()
	health := mat.NewVecDense(r, nil)
	health.MulVec(hf, activity)
	return health, nil
}

// HealthMatrix returns spatially- and industry-explicit air quality-related health impacts caused by the
// specified economic demand. In the result matrix, the rows represent air quality
// model grid cells and the columns represent emitters.
func (e *SpatialEIO) HealthMatrix(ctx context.Context, request *eieiorpc.HealthMatrixInput) (*eieiorpc.Matrix, error) {
	hr, ok := e.hr[request.HR]
	if !ok {
		return nil, fmt.Errorf("eieio: hazard ratio function `%s` is not registered", request.HR)
	}
	hf, err := e.healthFactors(ctx, request.AQM, Pollutant(request.Pollutant), request.Population, Year(request.Year), hr) // rows = grid cells, cols = industries
	if err != nil {
		return nil, err
	}

	activity, err := e.economicImpactsSCC(array2vec(request.Demand.Data), Year(request.Year), Location(request.Location)) // rows = industries
	if err != nil {
		return nil, err
	}

	r, c := hf.Dims()
	health := mat.NewDense(r, c, nil)
	health.Apply(func(_, j int, v float64) float64 {
		// Multiply each emissions factor column by the corresponding activity row.
		return v * activity.At(j, 0)
	}, hf)
	return mat2rpc(health), nil
}

type concPolPopYearHRAQM struct {
	pol  Pollutant
	year Year
	pop  string
	hr   epi.HRer
	aqm  string
}

// healthFactors returns spatially-explicit air quality-related health impacts per unit of economic
// production for each industry. In the result matrix, the rows represent
// air quality model grid cells and the columns represent industries.
func (e *SpatialEIO) healthFactors(ctx context.Context, aqm string, pol Pollutant, pop string, year Year, HR epi.HRer) (*mat.Dense, error) {
	e.loadHealthFactorsOnce.Do(func() {
		e.healthFactorCache = loadCacheOnce(e.healthFactorsWorker, 1, 1, e.EIEIOCache, matrixMarshal, matrixUnmarshal)
	})
	key := fmt.Sprintf("healthFactors_%s_%v_%v_%d_%s", aqm, pol, pop, year, HR.Name())
	rr := e.healthFactorCache.NewRequest(ctx, concPolPopYearHRAQM{pol: pol, year: year, pop: pop, hr: HR, aqm: aqm}, key)
	resultI, err := rr.Result()
	if err != nil {
		return nil, fmt.Errorf("bea: healthFactors: %s: %v", key, err)
	}
	return resultI.(*mat.Dense), nil
}

// healthFactorsWorker returns spatially-explicit pollution concentrations per unit of economic
// production for each industry. In the result matrix, the rows represent
// air quality model grid cells and the columns represent industries.
func (e *SpatialEIO) healthFactorsWorker(ctx context.Context, request interface{}) (interface{}, error) {
	polyearHR := request.(concPolPopYearHRAQM)
	prod, err := e.domesticProductionSCC(polyearHR.year)
	if err != nil {
		return nil, err
	}
	var healthFac *mat.Dense
	for i, refTemp := range e.SpatialRefs {
		if len(refTemp.SCCs) == 0 {
			return nil, fmt.Errorf("bea: industry %d; no SCCs", i)
		}
		ref := refTemp
		ref.EmisYear = int(polyearHR.year)
		ref.AQM = polyearHR.aqm
		health, err := e.CSTConfig.HealthSurrogate(ctx, &ref, polyearHR.hr.Name())
		if err != nil {
			return nil, err
		}
		var industryHealth *sparse.DenseArray
		var pol string
		switch polyearHR.pol {
		case PNH4:
			pol = "pNH4"
		case PNO3:
			pol = "pNO3"
		case PSO4:
			pol = "pSO4"
		case SOA:
			pol = "SOA"
		case PrimaryPM25:
			pol = "PrimaryPM25"
		case TotalPM25:
			pol = "TotalPM2_5"
		default:
			return nil, fmt.Errorf("eieio.health: invalid pollutant %v", polyearHR.pol)
		}
		if _, ok := health[polyearHR.pop]; !ok {
			return nil, fmt.Errorf("eieio.health: invalid population %v", polyearHR.pop)
		}
		industryHealth, ok := health[polyearHR.pop][pol]
		if !ok {
			return nil, fmt.Errorf("eieio.health: invalid pollutant %v", pol)
		}
		if i == 0 {
			healthFac = mat.NewDense(industryHealth.Shape[0], len(e.SpatialRefs), nil)
		}
		for r, v := range industryHealth.Elements {
			// The health factor is the industry health impacts divided by the
			// industry economic production.
			healthFac.Set(r, i, v/prod.At(i, 0))
		}
	}
	return healthFac, nil
}
