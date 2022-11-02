/*
Copyright © 2022 the InMAP-Canada authors.
This file is part of InMAP-Canada.

InMAP is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

InMAP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with InMAP.  If not, see <http://www.gnu.org/licenses/>.
*/

package inmap

import (
	"fmt"
	//"math"
	"time"

	"github.com/ctessum/atmos/seinfeld"
	"github.com/ctessum/atmos/wesely1989"

	"github.com/ctessum/sparse"
	//TR added 20220422 for aSOA/bSOA partition coefficient estimation. Could also copy to keep internal.
	//"github.com/sajari/regression"
)

//TR - update to GEMMACH variables
// GEM-MACH variables currently used:
/* hc5,hc8,olt,oli,tol,xyl,csl,cvasoa1,cvasoa2,cvasoa3,cvasoa4,iso,api,sesq,lim,
cvbsoa1,cvbsoa2,cvbsoa3,cvbsoa4,asoa1i,asoa1j,asoa2i,asoa2j,asoa3i,asoa3j,asoa4i,
asoa4j,bsoa1i,bsoa1j,bsoa2i,bsoa2j,bsoa3i,bsoa3j,bsoa4i,bsoa4j,no,no2,no3ai,no3aj,
so2,sulf,so4ai,so4aj,nh3,nh4ai,nh4aj,PM2_5_DRY,U,V,W,PBLH,PH,PHB,HFX,UST,PBLH,T,
PB,P,ho,h2o2,LU_INDEX,QRAIN,CLDFRA,QCLOUD,ALT,SWDOWN,GLW */

//const wrfFormat = "2006-01-02_15_04_05"
//20220907 Changed to GEMMACH time format. Pulled from test netcdf file provided by ECCC
const gemFormat = "2006-01-02T15:04:05.000000000"

// GEMMACH is an InMAP preprocessor for GEM-MACH output.
type GEMMACH struct {
	aVOC, bVOC, tSOA, aSOA, bSOA, nox, no, no2, pNO, sox, pS, nh3, pNH, totalPM25, ho, h2o2 map[string]float64

	start, end time.Time

	gemOut string

	recordDelta, fileDelta time.Duration

	msgChan chan string
}

// NewGEMMACH initializes a GEM-MACH preprocessor from the given
// configuration information.
// gemOut is the location of GEM-MACH output files.
// [DATE] should be used as a wild card for the simulation date.
// startDate and endDate are the dates of the beginning and end of the
// simulation, respectively, in the format "YYYYMMDD".
// If msgChan is not nil, status messages will be sent to it.
func NewGEMMACH(gemOut, startDate, endDate string, msgChan chan string) (*GEMMACH, error) {
	w := GEMMACH{
		// These maps contain the GEM-MACH variables that make
		// up the chemical species groups, as well as the
		// multiplication factors required to convert concentrations
		// to mass fractions [μg/kg dry air].

		// RACM VOC species and molecular weights (g/mol);
		// Only includes anthropogenic precursors to SOA from
		// anthropogenic (aSOA) and biogenic (bSOA) sources as
		// in Ahmadov et al. (2012)
		// Assume condensible vapor from SOA has molar mass of 70
		// TR - Need to adapt to GEMMACH aVOC variables - to be confirmed by ECCC
		// TR - Currently assuming TA3/TA2 are aVOCs
		//TRNotes - map is like a dict - maps keys to variables.
		//Right now, all we have is TTOL
		aVOC: map[string]float64{
			"TTOL": 1.,
		},
		//"TTOL": 1,"TA3": 1,"TA2": 1,

		//TR - Need to adapt to GEMMACH aVOC variables - to be confirmed by ECCC
		bVOC: map[string]float64{
			"TTOL": 0.5,
		},
		//"TISO": 1, "TPIN": 1, "TSES": 1,

		// For aSOA and bSOA, we will estimate the partition coefficient from a regression vs the aVOC/bVOC species
		//and the total PM2.5 SOA (TOC1)
		tSOA: map[string]float64{"TOC1": 1.},
		//May have to declare? Should pick up that they are arrays of floats. This may not work if all
		//values are read in one-by-one, rather than as vectors/arrays
		//Kpa,Kpb := Kpest(tSOA,aVOC,bVOC)

		//aSOA - calculate from TSOA
		//aSOA: map[string]float64{"aSOA": SOAest(tSOA, aVOC, bVOC, true)}, //Kpest(tSOA,aVOC,bVOC)[0]
		aSOA: map[string]float64{"TOC1": 0.5},
		// VBS SOA species (biogenic only) [μg/kg dry air].
		//bSOA: map[string]float64{"bSOA": SOAest(tSOA, aVOC, bVOC, false)},
		bSOA: map[string]float64{"TOC1": 0.5},
		// NOx is RACM NOx species. We are only interested in the mass
		// of Nitrogen, rather than the mass of the whole molecule, so
		// we use the molecular weight of Nitrogen.
		nox: map[string]float64{"TNO": 1., "TNO2": 1.},
		// pNO is the Nitrogen fraction of MADE particulate
		// NO species [μg/kg dry air].
		pNO: map[string]float64{"TNI1": 1., "TNO3": 1.},
		// SOx is the RACM SOx species. We are only interested in the mass
		// of Sulfur, rather than the mass of the whole molecule, so
		// we use the molecular weight of Sulfur.
		//TR/SB - converted from PPB to PPM by dividing by 1000
		sox: map[string]float64{"S2": ppmvToUgKg(mwS) / 1000.0},
		// pS is the Sulfur fraction of the MADE particulate
		// Sulfur species [μg/kg dry air].
		pS: map[string]float64{"TSU1": 1.},
		// NH3 is ammonia. We are only interested in the mass
		// of Nitrogen, rather than the mass of the whole molecule, so
		// we use the molecular weight of Nitrogen.
		nh3: map[string]float64{"TNH3": 1.},
		// pNH is the Nitrogen fraction of the MADE particulate
		// ammonia species [μg/m³].
		pNH: map[string]float64{"TAM1": 1.},
		// totalPM25 is total mass of PM2.5  [μg/m3].
		totalPM25: map[string]float64{"AF": 1.},
		// Hydroxy radical is concentratoon of the hydroxy radical.
		//Convert from ug/kg to ppmV for consistency
		ho: map[string]float64{"TOH": UgKgToppmv(mwOH)},
		//hydrogen peroxide is concentratoon of H202.
		//Convert from ug/kg to ppmV for consistency
		h2o2: map[string]float64{"TH22": UgKgToppmv(mwH2O2)},

		gemOut:  gemOut,
		msgChan: msgChan,
	}

	var err error
	w.start, err = time.Parse(inDateFormat, startDate)
	if err != nil {
		return nil, fmt.Errorf("inmap: GEM-MACH preprocessor start time: %v", err)
	}
	w.end, err = time.Parse(inDateFormat, endDate)
	if err != nil {
		return nil, fmt.Errorf("inmap: GEM-MACH preprocessor end time: %v", err)
	}

	if !w.end.After(w.start) {
		if err != nil {
			return nil, fmt.Errorf("inmap: GEM-MACH preprocessor end time %v is not after start time %v", w.end, w.start)
		}
	}

	w.recordDelta, err = time.ParseDuration("1h")
	if err != nil {
		return nil, fmt.Errorf("inmap: GEM-MACH preprocessor recordDelta: %v", err)
	}
	w.fileDelta, err = time.ParseDuration("24h")
	if err != nil {
		return nil, fmt.Errorf("inmap: GEM-MACH preprocessor fileDelta: %v", err)
	}
	return &w, nil
}

// ppmvToUgKg returns a multiplier to convert a concentration in
// ppmv dry air to a mass fraction [micrograms per kilogram dry air]
// for a chemical species with the given molecular weight in g/mol.
//func ppmvToUgKg(mw float64) float64 {
//return mw * 1000.0 / MWa
//}
func UgKgToppmv(mw float64) float64 {
	return mw * MWa / 1000.0
}

//TR added for linear regression of aSOA/bSOA
//This is more like pseudo-code, need to figure out how to get all time
//values at each spatial grid cell for the regression
//Need to add a for loop to go through each of the arrays in time
/*
func SOAest(tSOA []float64, aVOC, bVOC float64, return_aSOA bool) []float64 {
	//var mat [][]float64
	r := new(regression.Regression)
	r.SetObserved("tSOA")
	r.SetVar(0, "aVOC")
	r.SetVar(1, "bVOC")
	//Everything has to have the same length.
	//dpts := append(tSOA, aVOC, bVOC)
	veclen := len(tSOA)
	mat := [veclen][3]float64{tSOA, aVOC, bVOC}
	dpts := regression.MakeDataPoints(mat, 0)
	//var dpts []float64
	//for
	r.Train(dpts)
	r.Run()
	Kpa := r.Coeff(0)
	Kpb := r.Coeff(1)
	if return_aSOA == true {
		SOA := Kpa * aVOC
	} else {
		SOA := Kpb * bVOC
	}
	//aSOA := Kpa*aVOC
	//bSOA := Kpb*bVOC
	return SOA
}
*/
/*
//The ellipses allows for an arbitrary number of variables (vars)

func linest(obs float64, vars ...float64) float64 {
	r := new(regression.Regression)
	ind := 0
	for j := range vars {
        r.SetVar(ind, "Inhabitants")
		ind += 1
    }
	c1 := r.Coeff(0)
	c2 := r.Coeff(1)
	return mw * 1000.0 / MWa
}
*/
//Syntax for these methods since I find them a bit confusing. Basically, a method is a function for a struct (or similar)
//These functions set up w as a "receiver" of type GEMMACH (*GEMMACH sets up a pointer to the GEMMACH struct)that then has the method  "read" or whatever.
//So it goes func (receiver *struct) funcname(input inputtype) NextData (which is a function to get the next timestep data from preproc.go) {stuff the function does}.
//This format lets you call the function as w.read(varName) etc.

func (w *GEMMACH) read(varName string) NextData {
	return nextDataNCF(w.gemOut, wrfFormat, varName, w.start, w.end, w.recordDelta, w.fileDelta, readNCF, w.msgChan)
}

func (w *GEMMACH) readGroupAlt(varGroup map[string]float64) NextData {
	return nextDataGroupAltNCF(w.gemOut, gemFormat, varGroup, w.ALT(), w.start, w.end, w.recordDelta, w.fileDelta, readNCF, w.msgChan)
}

func (w *GEMMACH) readGroup(varGroup map[string]float64) NextData {
	return nextDataGroupNCF(w.gemOut, gemFormat, varGroup, w.start, w.end, w.recordDelta, w.fileDelta, readNCF, w.msgChan)
}

//For these three functions - changed from ALT to AF as GEMMACH
//outputs do not contain ALT, but do contain AF
// Nx helps fulfill the Preprocessor interface by returning
// the number of grid cells in the West-East direction.
func (w *GEMMACH) Nx() (int, error) {
	f, ff, err := ncfFromTemplate(w.gemOut, gemFormat, w.start)
	if err != nil {
		return -1, fmt.Errorf("nx: %v", err)
	}
	defer f.Close()
	return ff.Header.Lengths("AF")[3], nil
}

// Ny helps fulfill the Preprocessor interface by returning
// the number of grid cells in the South-North direction.
func (w *GEMMACH) Ny() (int, error) {
	f, ff, err := ncfFromTemplate(w.gemOut, gemFormat, w.start)
	if err != nil {
		return -1, fmt.Errorf("ny: %v", err)
	}
	defer f.Close()
	return ff.Header.Lengths("AF")[2], nil
}

// Nz helps fulfill the Preprocessor interface by returning
// the number of grid cells in the below-above direction.
func (w *GEMMACH) Nz() (int, error) {
	f, ff, err := ncfFromTemplate(w.gemOut, gemFormat, w.start)
	if err != nil {
		return -1, fmt.Errorf("nz: %v", err)
	}
	defer f.Close()
	return ff.Header.Lengths("AF")[1], nil
}

// PBLH helps fulfill the Preprocessor interface by returning
// planetary boundary layer height [m].
func (w *GEMMACH) PBLH() NextData { return w.read("H") }

// Height helps fulfill the Preprocessor interface by returning
// layer heights above ground level calculated based on geopotential height.
// For more information, refer to
// http://www.openwfm.org/wiki/How_to_interpret_WRF_variables.
//TRSB - Used geopotential height directly for height, converted from decametres
func (w *GEMMACH) Height() NextData { return w.read("GZ") }

//func (w *GEMMACH) Height() NextData { return w.read("GZ") *10.0 }

/*TRSB
func (w *GEMMACH) Height() NextData {
	// ph is perturbation geopotential height [m/s2].
	phFunc := w.read("PH")
	// phb is baseline geopotential height [m/s2].
	phbFunc := w.read("PHB")
	return func() (*sparse.DenseArray, error) {
		ph, err := phFunc()
		if err != nil {
			return nil, err
		}
		phb, err := phbFunc()
		if err != nil {
			return nil, err
		}
		return geopotentialToHeight(ph, phb), nil
	}
}

func geopotentialToHeight(ph, phb *sparse.DenseArray) *sparse.DenseArray {
	layerHeights := sparse.ZerosDense(ph.Shape...)
	for k := 0; k < ph.Shape[0]; k++ {
		for j := 0; j < ph.Shape[1]; j++ {
			for i := 0; i < ph.Shape[2]; i++ {
				h := (ph.Get(k, j, i) + phb.Get(k, j, i) -
					ph.Get(0, j, i) - phb.Get(0, j, i)) / g // m
				layerHeights.Set(h, k, j, i)
			}
		}
	}
	return layerHeights
}
*/
// ALT helps fulfill the Preprocessor interface by returning
// inverse air density [m3/kg].
func (w *GEMMACH) ALT() NextData {
	rhoFunc := w.read("RHO") // Density of air kg m-2
	return func() (*sparse.DenseArray, error) {
		rho, err := rhoFunc()
		if err != nil {
			return nil, err
		}
		alt := sparse.ZerosDense(rho.Shape...)
		for i, rhoV := range rho.Elements {
			// molec HO / cm3 * m3 / kg air * kg air/molec. air* cm3/m3 * ppm
			alt.Elements[i] = 1 / rhoV
		}
		//mult := sparse.ZerosDense(rho.Shape...)
		//alt := mult / rho
		return alt, nil
	}
}

// U helps fulfill the Preprocessor interface by returning
//Converting from knots to m/s
// West-East wind speed [m/s].
func (w *GEMMACH) U() NextData {
	uufunc := w.read("UU")
	return func() (*sparse.DenseArray, error) {
		UU, err := uufunc()
		if err != nil {
			return nil, err
		}
		//UU = UU.Scale(463.0 / 900.0)
		return UU.ScaleCopy(463.0 / 900.0), nil
	}
}

// V helps fulfill the Preprocessor interface by returning
// South-North wind speed [m/s].
//func (w *GEMMACH) V() NextData { return w.read("VV") * 463.0 / 900.0 }
func (w *GEMMACH) V() NextData {
	vvfunc := w.read("VV")
	return func() (*sparse.DenseArray, error) {
		VV, err := vvfunc()
		if err != nil {
			return nil, err
		}
		//UU = UU.Scale(463.0 / 900.0)
		return VV.ScaleCopy(463.0 / 900.0), nil
	}
}

// W helps fulfill the Preprocessor interface by returning
// below-above wind speed [m/s].
/*
func (w *GEMMACH) W() NextData {
	WW := w.read("WW")
	Height := w.Height()
	P := w.P()
	return convert_vertspeed(WW, Height, P)
}
//Convert vertical windspeed from Pa/s to m/s
func convert_vertspeed(WW, GZ, P *sparse.DenseArray) *sparse.DenseArray {
	vertspeeds := sparse.ZerosDense(WW.Shape...)
	for k := 0; k < WW.Shape[0]; k++ {
		for j := 0; j < WW.Shape[1]; j++ {
			for i := 0; i < WW.Shape[2]; i++ {
				slope := (GZ.Get(k+1, j, i) - GZ.Get(k, j, i)) /
					(P.Get(k+1, j, i) - P.Get(k, j, i))
				WWprime := WW.Get(k, j, i) * slope
				vertspeeds.Set(WWprime, k, j, i)
			}
		}
	}
	return vertspeeds
}
*/
func (w *GEMMACH) W() NextData {
	wwFunc := w.read("WW") // windspeed  in pa/s
	gzFunc := w.Height()   //Geopotential height (m)
	ppFunc := w.P()
	return func() (*sparse.DenseArray, error) {
		WW, err := wwFunc()
		if err != nil {
			return nil, err
		}
		GZ, err := gzFunc()
		if err != nil {
			return nil, err
		}
		PP, err := ppFunc()
		if err != nil {
			return nil, err
		}
		out := sparse.ZerosDense(WW.Shape...)
		for k := 0; k < WW.Shape[0]; k++ {
			for j := 0; j < WW.Shape[1]; j++ {
				for i := 0; i < WW.Shape[2]; i++ {
					slope := (GZ.Get(k+1, j, i) - GZ.Get(k, j, i)) /
						(PP.Get(k+1, j, i) - PP.Get(k, j, i))
					WWprime := WW.Get(k, j, i) * slope
					out.Set(WWprime, k, j, i)
				}
			}
		}
		return out, nil
	}
}

// AVOC helps fulfill the Preprocessor interface.
func (w *GEMMACH) AVOC() NextData { return w.readGroupAlt(w.aVOC) }

// BVOC helps fulfill the Preprocessor interface.
func (w *GEMMACH) BVOC() NextData { return w.readGroupAlt(w.bVOC) }

// NOx helps fulfill the Preprocessor interface.
func (w *GEMMACH) NOx() NextData { return w.readGroupAlt(w.nox) }

// SOx helps fulfill the Preprocessor interface.
func (w *GEMMACH) SOx() NextData { return w.readGroupAlt(w.sox) }

// NH3 helps fulfill the Preprocessor interface.
func (w *GEMMACH) NH3() NextData { return w.readGroupAlt(w.nh3) }

// ASOA helps fulfill the Preprocessor interface.
func (w *GEMMACH) ASOA() NextData { return w.readGroupAlt(w.aSOA) }

// BSOA helps fulfill the Preprocessor interface.
func (w *GEMMACH) BSOA() NextData { return w.readGroupAlt(w.bSOA) }

// PNO helps fulfill the Preprocessor interface.
func (w *GEMMACH) PNO() NextData { return w.readGroupAlt(w.pNO) }

// PS helps fulfill the Preprocessor interface.
func (w *GEMMACH) PS() NextData { return w.readGroupAlt(w.pS) }

// PNH helps fulfill the Preprocessor interface.
//TRSB - changed to "read group" not "read group alt" as value already in ug/kg
func (w *GEMMACH) PNH() NextData { return w.readGroup(w.pNH) }

// TotalPM25 helps fulfill the Preprocessor interface.
func (w *GEMMACH) TotalPM25() NextData { return w.readGroup(w.totalPM25) }

// SurfaceHeatFlux helps fulfill the Preprocessor interface
// by returning heat flux at the surface [W/m2].
//FC5 is the aggregate heat flux
func (w *GEMMACH) SurfaceHeatFlux() NextData { return w.read("FC5") }

// UStar helps fulfill the Preprocessor interface
// by returning friction velocity [m/s].
func (w *GEMMACH) UStar() NextData { return w.read("UE") }

// T helps fulfill the Preprocessor interface by
// returning temperature [K].
//func (w *GEMMACH) T() NextData { return 273.15 + w.read("TT") }

func (w *GEMMACH) T() NextData {
	TTFunc := w.read("TT") // Density of air kg m-2
	return func() (*sparse.DenseArray, error) {
		TT, err := TTFunc()
		if err != nil {
			return nil, err
		}
		//alt := sparse.ZerosDense(TT.Shape...)
		for i, TTV := range TT.Elements {
			// molec HO / cm3 * m3 / kg air * kg air/molec. air* cm3/m3 * ppm
			TT.Elements[i] = TTV + 273.15
		}
		//mult := sparse.ZerosDense(rho.Shape...)
		//alt := mult / rho
		return TT, nil
	}
}

/*
thetaFunc := w.read("T") // perturbation potential temperature [K]
	pFunc := w.P()           // Pressure [Pa]
	return wrfTemperatureConvert(thetaFunc, pFunc)
}

func wrfTemperatureConvert(thetaFunc, pFunc NextData) NextData {
	return func() (*sparse.DenseArray, error) {
		thetaPerturb, err := thetaFunc() // perturbation potential temperature [K]
		if err != nil {
			return nil, err
		}
		p, err := pFunc() // Pressure [Pa]
		if err != nil {
			return nil, err
		}

		T := sparse.ZerosDense(thetaPerturb.Shape...)
		for i, tp := range thetaPerturb.Elements {
			T.Elements[i] = thetaPerturbToTemperature(tp, p.Elements[i])
		}
		return T, nil
	}
}

// thetaPerturbToTemperature converts perburbation potential temperature
// to ambient temperature for the given pressure (p [Pa]).
func thetaPerturbToTemperature(thetaPerturb, p float64) float64 {
	const (
		po    = 101300. // Pa, reference pressure
		kappa = 0.2854  // related to von karman's constant
	)
	pressureCorrection := math.Pow(p/po, kappa)
	// potential temperature, K
	θ := thetaPerturb + 300.
	// Ambient temperature, K
	return θ * pressureCorrection
}
*/
// P helps fulfill the Preprocessor interface
// by returning pressure [Pa] as pressure level converted
// from atmosphere to Pa
//func (w *GEMMACH) P() NextData { return 133125.0 * w.read("level1") }
func (w *GEMMACH) P() NextData {
	ppfunc := w.read("level1")
	return func() (*sparse.DenseArray, error) {
		PP, err := ppfunc()
		if err != nil {
			return nil, err
		}
		//UU = UU.Scale(463.0 / 900.0)
		return PP.ScaleCopy(133125.0), nil
	}
}

/*
func (w *GEMMACH) P() NextData {
	pbFunc := w.read("PB") // baseline pressure [Pa]
	pFunc := w.read("P")   // perturbation pressure [Pa]
	return wrfPressureConvert(pFunc, pbFunc)
}

func wrfPressureConvert(pFunc, pbFunc NextData) NextData {
	return func() (*sparse.DenseArray, error) {
		pb, err := pbFunc() // baseline pressure [Pa]
		if err != nil {
			return nil, err
		}
		p, err := pFunc() // perturbation pressure [Pa]
		if err != nil {
			return nil, err
		}
		P := pb.Copy()
		P.AddDense(p)
		return P, nil
	}
}
*/
// HO helps fulfill the Preprocessor interface
// by returning hydroxyl radical concentration [ppmv].
func (w *GEMMACH) HO() NextData { return w.readGroup(w.ho) }

// H2O2 helps fulfill the Preprocessor interface
// by returning hydrogen peroxide concentration [ppmv].
func (w *GEMMACH) H2O2() NextData { return w.readGroup(w.h2o2) }

// SeinfeldLandUse helps fulfill the Preprocessor interface
// by returning land use categories as
// specified in github.com/ctessum/atmos/seinfeld.
func (w *GEMMACH) SeinfeldLandUse() NextData {
	luFunc := w.read("LU_INDEX") // USGS land use index
	return GEMSeinfeldLandUse(luFunc)
}

func GEMSeinfeldLandUse(luFunc NextData) NextData {
	return func() (*sparse.DenseArray, error) {
		lu, err := luFunc() // USGS land use index
		if err != nil {
			return nil, err
		}
		o := sparse.ZerosDense(lu.Shape...)
		for j := 0; j < lu.Shape[0]; j++ {
			for i := 0; i < lu.Shape[1]; i++ {
				o.Set(float64(GEMseinfeld[f2i(lu.Get(j, i))]), j, i)
			}
		}
		return o, nil
	}
}

// ***Need to convery vegetation fraction to land use code - take whatever land use is max and call it VF_INDEX.***
// GEMseinfeld lookup table to go from GEM land classes to land classes for
// particle dry deposition.
//TR - commented land use is WRF land use definition, left hand side is INMAP definition.
var GEMseinfeld = []seinfeld.LandUseCategory{
	seinfeld.Shrubs,    // Mixed shrubs VF01
	seinfeld.Deciduous, //'Mixed Forest' Mixed wood forest VF02
	seinfeld.Desert,    //'Barren or Sparsely Vegetated' desert VF03
	seinfeld.Deciduous, //'Wooded Wetland' Swamp VF04
	seinfeld.Desert,    //'Barren or Sparsely Vegetated' Tundra VF05
	seinfeld.Desert,    //'Urban and Built-Up Land' Urban VF06
	seinfeld.Grass,     //'Irrigated Cropland and Pasture' Irrigated crops VF07
	seinfeld.Grass,     // Cotton VF08
	seinfeld.Grass,     // Maize VF09
	seinfeld.Grass,     // Sugar VF10
	seinfeld.Grass,     // Rice VF11
	seinfeld.Grass,     // Crops VF12
	seinfeld.Grass,     // Long Grass VF13
	seinfeld.Shrubs,    // 'Mixed Shrubland/Grassland' Short grass and forbs VF14
	seinfeld.Shrubs,    // Thorn shrubs VF15
	seinfeld.Shrubs,    // Deciduous shrubs VF16
	seinfeld.Shrubs,    // evergreen broadleaf shrubs VF17
	seinfeld.Deciduous, // drought deciduous trees VF18
	seinfeld.Deciduous, // Tropical broadleaf trees VF19
	seinfeld.Deciduous, // deciduous broadleaf trees VF20
	seinfeld.Evergreen, // deciduous needle-leaf trees VF21
	seinfeld.Deciduous, // Evergreen broadleaf trees VF22
	seinfeld.Evergreen, // Evergreen needle-leaf trees VF23
	seinfeld.Desert,    // Inland lake VF24
	seinfeld.Desert,    // Ice VF25
	seinfeld.Desert,    // Water VF26
}

// WeselyLandUse helps fulfill the Preprocessor interface
// by returning land use categories as
// specified in github.com/ctessum/atmos/wesely1989.
func (w *GEMMACH) WeselyLandUse() NextData {
	luFunc := w.read("VF_INDEX") // USGS land use index
	return GEMWeselyLandUse(luFunc)
}

func GEMWeselyLandUse(luFunc NextData) NextData {
	return func() (*sparse.DenseArray, error) {
		lu, err := luFunc() // USGS land use index
		if err != nil {
			return nil, err
		}
		o := sparse.ZerosDense(lu.Shape...)
		for j := 0; j < lu.Shape[0]; j++ {
			for i := 0; i < lu.Shape[1]; i++ {
				o.Set(float64(GEMVFwesely[f2i(lu.Get(j, i))]), j, i)
			}
		}
		return o, nil
	}
}

// GEMVFwesely lookup table to go from GEM vegetation fraction to land classes for
// gas dry deposition.
//TR - need to convert from GEMMACH wesely to INMAP wesely
var GEMVFwesely = []wesely1989.LandUseCategory{
	wesely1989.RangeAg,     //'Mixed Shrubland/Grassland' mixed_shrubs
	wesely1989.MixedForest, //'Mixed Forest' mixed wood forest VF02
	wesely1989.Barren,      //'White Sand' Desert VF03
	wesely1989.Wetland,     //'Wooded Wetland' Swamp VF04
	wesely1989.RockyShrubs, //'Mixed Tundra' Tundra VF05
	wesely1989.Urban,       //'Urban and Built-Up Land' Urban VF06
	wesely1989.RangeAg,     //'Irrigated Cropland and Pasture' irrigated crops VF07
	wesely1989.RangeAg,     //'Irrigated Cropland and Pasture' Cotton VF08
	wesely1989.RangeAg,     //'Dryland Cropland and Pasture' Maize VF09
	wesely1989.RangeAg,     //'Irrigated Cropland and Pasture' Sugar VF10
	wesely1989.RangeAg,     //'Irrigated Cropland and Pasture' Rice VF11
	wesely1989.RangeAg,     //'Cropland/Grassland Mosaic' Crops VF12
	wesely1989.Range,       //'Grassland' Long grass VF13
	wesely1989.RangeAg,     //'Mixed Shrubland/Grassland' short grass/forbs VF14
	wesely1989.RockyShrubs, //'Shrubland' Thorn shrubs VF15
	wesely1989.RockyShrubs, //'Shrubland' Deciduous shrubs VF16
	wesely1989.RockyShrubs, //'Shrubland'  Evergreen broadleaf shrubs VF 17
	wesely1989.Deciduous,   //'Deciduous Broadleaf Forest' drought deciduous trees VF18
	wesely1989.Deciduous,   //'Evergreen Broadleaf Forest' Tropical broadleaf trees VF19
	wesely1989.Deciduous,   //'Deciduous Broadleaf Forest' deciduous broadleaf tree VF20
	wesely1989.Coniferous,  //'Deciduous Needleleaf Forest' deciduous needleleaf trees VF21
	wesely1989.Deciduous,   //'Evergreen Broadleaf Forest' Evergreen broadleaf trees VF22
	wesely1989.Coniferous,  //'Evergreen Needleleaf Forest' evergreen needleleaf trees VF23
	wesely1989.Water,       //'Water Bodies' inland lake VF24
	wesely1989.Barren,      //'Snow or Ice' ice VF25
	wesely1989.Water,       //'Water Bodies' water VF26
}

// Z0 helps fulfill the Preprocessor interface by
// returning roughness length.
func (w *GEMMACH) Z0() NextData {
	LUIndexFunc := w.read("VF_INDEX") //GEM VF land use index
	return GEMZ0(LUIndexFunc)
}

func GEMZ0(LUIndexFunc NextData) NextData {
	return func() (*sparse.DenseArray, error) {
		luIndex, err := LUIndexFunc()
		if err != nil {
			return nil, err
		}
		zo := sparse.ZerosDense(luIndex.Shape...)
		for i, lu := range luIndex.Elements {
			zo.Elements[i] = GEMz0[f2i(lu)] // roughness length [m]
		}
		return zo, nil
	}
}

// GEMz0 holds Roughness lengths for GEM land classes ([m]), from WRF file
// TR - find GEMMACH roughness length definitions.
var GEMz0 = []float64{0.001, 0.0003, 0.001, 1.5, 3.5, 1,
	2.0, 3.0, .80, .05, .15, .15, .02, .08, 0.08, .08, .35,
	.25, 0.1, .08, 1.35, .01, .05, .05, 1.5, .05}

// QRain helps fulfill the Preprocessor interface by
// returning rain mass fraction.
//TR -Sahil to do
//func (w *GEMMACH) QRain() NextData { return w.read("QRAIN") }
// QRain helps fulfill the Preprocessor interface by returning
// rain mass fraction based on the GEM precipitation rate [kg m-2 s-1]
// and the assumption (from the EMEP model wet deposition algorithm)
// that raindrops are falling at 5 m/s.
func (w *GEMMACH) QRain() NextData {
	//PFLCUFunc := w.read("PR")/8760.0     // Quantity of precipitation (m/hr)
	//PFLLSanFunc := w.read("PC")/8760.0 // Implicit accumulated precipitation (m/hr)
	PrecipFunc := w.read("RT") // Total precipitation rate (m/s)
	altFunc := w.ALT()         // 1/density
	//PFLLSanFunc := w.read("VDR")       // Deposition velocity (m/s) - not possible to find
	return func() (*sparse.DenseArray, error) {
		/*
			pflcu, err := PFLCUFunc()
			if err != nil {
				return nil, err
			}
			pfllSan, err := PFLLSanFunc()
			if err != nil {
				return nil, err
			}
			}
		*/

		precipgem, err := PrecipFunc()
		if err != nil {
			return nil, err
		}
		alt, err := altFunc()
		if err != nil {
			return nil, err
		}
		//}
		const Vdr = 5.0 // droplet velocity [m/s] EMEP wet dep algorithm assumption
		qRain := sparse.ZerosDense(alt.Shape...)
		// pflcu and pfllSan are staggered but qRain is not, so
		// we average the rain flux values (pflcu and pfllSan) at the
		// bottom and top of each grid cell.
		for k := 0; k < qRain.Shape[0]; k++ {
			for j := 0; j < qRain.Shape[1]; j++ {
				for i := 0; i < qRain.Shape[2]; i++ {
					// From EMEP algorithm: P = QRAIN * Vdr * ρgas => QRAIN = P / Vdr / ρgas
					// [kg m-2 s-1] / [m s-1] * [m3 kg-1]
					if k == 0 {
						//q := (pflcu.Get(k, j, i) + pfllSan.Get(k, j, i)) / Vdr
						q := (precipgem.Get(k, j, i)) / Vdr
						qRain.Set(q, k, j, i)
					} else {
						q := 0.0
						qRain.Set(q, k, j, i)
					}
				}
			}
		}
		return qRain, nil
	}
}

// CloudFrac helps fulfill the Preprocessor interface
// by returning the fraction of each grid cell filled
// with clouds [volume/volume].
func (w *GEMMACH) CloudFrac() NextData { return w.read("FN") }

// QCloud helps fulfill the Preprocessor interface by returning
// the mass fraction of cloud water in each grid cell [mass/mass].
func (w *GEMMACH) QCloud() NextData { return w.read("QC") }

// RadiationDown helps fulfill the Preprocessor interface by returning
// total downwelling radiation at ground level [W/m2].
func (w *GEMMACH) RadiationDown() NextData { return w.read("FB") }

/*	swDownFunc := w.read("SWDOWN") // downwelling short wave radiation at ground level [W/m2]
	glwFunc := w.read("GLW")       // downwelling long wave radiation at ground level [W/m2]
	return wrfRadiationDown(swDownFunc, glwFunc)
}

func wrfRadiationDown(swDownFunc, glwFunc NextData) NextData {
	return func() (*sparse.DenseArray, error) {
		swDown, err := swDownFunc() // downwelling short wave radiation at ground level [W/m2]
		if err != nil {
			return nil, err
		}
		glw, err := glwFunc() // downwelling long wave radiation at ground level [W/m2]
		if err != nil {
			return nil, err
		}
		rad := swDown.Copy()
		rad.AddDense(glw)
		return rad, nil
	}
}
*/
// SWDown helps fulfill the Preprocessor interface by returning
// downwelling short wave radiation at ground level [W/m2].
func (w *GEMMACH) SWDown() NextData { return w.read("FI") }

// GLW helps fulfill the Preprocessor interface by returning
// downwelling long wave radiation at                                                                                                                 ground level [W/m2].
//func (w *GEMMACH) GLW() NextData { return w.RadiationDown() - w.SWDown() }
func (w *GEMMACH) GLW() NextData {
	rdFunc := w.RadiationDown() // Density of air kg m-2
	swFunc := w.SWDown()
	return func() (*sparse.DenseArray, error) {
		RD, err := rdFunc()
		if err != nil {
			return nil, err
		}
		SW, err := swFunc()
		if err != nil {
			return nil, err
		}
		glw := sparse.ZerosDense(SW.Shape...)
		for i := range glw.Elements {
			// molec HO / cm3 * m3 / kg air * kg air/molec. air* cm3/m3 * ppm
			glw.Elements[i] = RD.Elements[i] + SW.Elements[i]
		}
		//mult := sparse.ZerosDense(rho.Shape...)
		//alt := mult / rho
		return glw, nil
	}
}
