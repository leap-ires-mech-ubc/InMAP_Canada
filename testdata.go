/*
Copyright © 2013 the InMAP authors.
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
along with InMAP.  If not, see <http://www.gnu.org/licenses/>.
*/

package inmap

import (
	"os"
	"strings"

	"github.com/ctessum/geom"
	"github.com/ctessum/geom/encoding/shp"
	"github.com/ctessum/sparse"
)

const (
	// TestGridSR is a spatial reference for testing.
	TestGridSR = `PROJCS["Lambert_Conformal_Conic_2SP",GEOGCS["GCS_unnamed ellipse",DATUM["D_unknown",SPHEROID["Unknown",6370997,0]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Lambert_Conformal_Conic_2SP"],PARAMETER["standard_parallel_1",33],PARAMETER["standard_parallel_2",45],PARAMETER["latitude_of_origin",40],PARAMETER["central_meridian",-97],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["Meter",1]]`
	// TestPopulationShapefile is a population shapefile for testing.
	TestPopulationShapefile = "tempPopulation.shp"
	// TestMortalityShapefile is mortality rate shapefile for testing.
	TestMortalityShapefile = "tempMortality.shp"
	// TestCTMDataFile is a CTMData file for testing.
	TestCTMDataFile = "tempCTMData.ncf"
)

// WriteTestPopShapefile writes out a population shapefile for testing.
var WriteTestPopShapefile = func() {
	// holder for test population data.
	type pop struct {
		geom.Polygon
		TotalPop, WhiteNoLat, Black, Native, Asian, Latino float64
	}

	// write out test population data.
	popData := []pop{
		// {
		// 	Polygon: []geom.Path{{
		// 		geom.Point{X: -3999, Y: -3999},
		// 		geom.Point{X: -3899, Y: -3999},
		// 		geom.Point{X: -3899, Y: -3899},
		// 		geom.Point{X: -3999, Y: -3899},
		// 		geom.Point{X: -3999, Y: -3999},
		// 	}},
		// 	TotalPop:   100000., // enough to split grid cell the smallest level
		// 	WhiteNoLat: 50000.,
		// 	Black:      20000.,
		// 	Asian:      8000.,
		// 	Native:     2000.,
		// 	Latino:     20000.,
		// },
		{
			Polygon: []geom.Path{{
				geom.Point{X: -31.815, Y: -39.4932},
				geom.Point{X: -18.765, Y: -39.4932},
				geom.Point{X: -18.765, Y: -23.1112},
				geom.Point{X: -31.815, Y: -23.1112},
				geom.Point{X: -31.815, Y: -39.4932},
			}},
			TotalPop:   100000., // enough to split grid cell the smallest level
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -31.815, Y: -23.1112},
				geom.Point{X: -18.765, Y: -23.1112},
				geom.Point{X: -18.765, Y: -6.72923},
				geom.Point{X: -31.815, Y: -6.72923},
				geom.Point{X: -31.815, Y: -23.1112},
			}},
			TotalPop:   100000., // enough to split grid cell the smallest level
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -31.815, Y: -6.72923},
				geom.Point{X: -18.765, Y: -6.72923},
				geom.Point{X: -18.765, Y: 9.65277},
				geom.Point{X: -31.815, Y: 9.65277},
				geom.Point{X: -31.815, Y: -6.72923},
			}},
			TotalPop:   100000., // enough to split grid cell the smallest level
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -31.815, Y: 9.65277},
				geom.Point{X: -18.765, Y: 9.65277},
				geom.Point{X: -18.765, Y: 26.03477},
				geom.Point{X: -31.815, Y: 26.03477},
				geom.Point{X: -31.815, Y: 9.65277},
			}},
			TotalPop:   100000., // enough to split grid cell the smallest level
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
		//Row 2
		{
			Polygon: []geom.Path{{
				geom.Point{X: -18.765, Y: -39.4932},
				geom.Point{X: -5.715, Y: -39.4932},
				geom.Point{X: -5.715, Y: -23.1112},
				geom.Point{X: -18.765, Y: -23.1112},
				geom.Point{X: -18.765, Y: -39.4932},
			}},
			TotalPop:   0., // Make sure not to split
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -18.765, Y: -23.1112},
				geom.Point{X: -5.715, Y: -23.1112},
				geom.Point{X: -5.715, Y: -6.72923},
				geom.Point{X: -18.765, Y: -6.72923},
				geom.Point{X: -18.765, Y: -23.1112},
			}},
			TotalPop:   0., // Make sure not to split
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -18.765, Y: -6.72923},
				geom.Point{X: -5.715, Y: -6.72923},
				geom.Point{X: -5.715, Y: 9.65277},
				geom.Point{X: -18.765, Y: 9.65277},
				geom.Point{X: -18.765, Y: -6.72923},
			}},
			TotalPop:   0., // Make sure not to split
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -18.765, Y: 9.65277},
				geom.Point{X: -5.715, Y: 9.65277},
				geom.Point{X: -5.715, Y: 26.03477},
				geom.Point{X: -18.765, Y: 26.03477},
				geom.Point{X: -18.765, Y: 9.65277},
			}},
			TotalPop:   0., // Make sure not to split
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
		//Row 3
		{
			Polygon: []geom.Path{{
				geom.Point{X: -5.715, Y: -39.4932},
				geom.Point{X: 7.335, Y: -39.4932},
				geom.Point{X: 7.335, Y: -23.1112},
				geom.Point{X: -5.715, Y: -23.1112},
				geom.Point{X: -5.715, Y: -39.4932},
			}},
			TotalPop:   0., // Make sure not to split
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -5.715, Y: -23.1112},
				geom.Point{X: 7.335, Y: -23.1112},
				geom.Point{X: 7.335, Y: -6.72923},
				geom.Point{X: -5.715, Y: -6.72923},
				geom.Point{X: -5.715, Y: -23.1112},
			}},
			TotalPop:   0., // Make sure not to split
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -5.715, Y: -6.72923},
				geom.Point{X: 7.335, Y: -6.72923},
				geom.Point{X: 7.335, Y: 9.65277},
				geom.Point{X: -5.715, Y: 9.65277},
				geom.Point{X: -5.715, Y: -6.72923},
			}},
			TotalPop:   0., // Make sure not to split
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -5.715, Y: 9.65277},
				geom.Point{X: 7.335, Y: 9.65277},
				geom.Point{X: 7.335, Y: 26.03477},
				geom.Point{X: -5.715, Y: 26.03477},
				geom.Point{X: -5.715, Y: 9.65277},
			}},
			TotalPop:   0., // Make sure not to split
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
		//Row 4
		{
			Polygon: []geom.Path{{
				geom.Point{X: 7.335, Y: -39.4932},
				geom.Point{X: 20.385, Y: -39.4932},
				geom.Point{X: 20.385, Y: -23.1112},
				geom.Point{X: 7.335, Y: -23.1112},
				geom.Point{X: 7.335, Y: -39.4932},
			}},
			TotalPop:   0., // Make sure not to split
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: 7.335, Y: -23.1112},
				geom.Point{X: 20.385, Y: -23.1112},
				geom.Point{X: 20.385, Y: -6.72923},
				geom.Point{X: 7.335, Y: -6.72923},
				geom.Point{X: 7.335, Y: -23.1112},
			}},
			TotalPop:   0., // Make sure not to split
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: 7.335, Y: -6.72923},
				geom.Point{X: 20.385, Y: -6.72923},
				geom.Point{X: 20.385, Y: 9.65277},
				geom.Point{X: 7.335, Y: 9.65277},
				geom.Point{X: 7.335, Y: -6.72923},
			}},
			TotalPop:   0., // Make sure not to split
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: 7.335, Y: 9.65277},
				geom.Point{X: 20.385, Y: 9.65277},
				geom.Point{X: 20.385, Y: 26.03477},
				geom.Point{X: 7.335, Y: 26.03477},
				geom.Point{X: 7.335, Y: 9.65277},
			}},
			TotalPop:   0., // Make sure not to split
			WhiteNoLat: 50000.,
			Black:      20000.,
			Asian:      8000.,
			Native:     2000.,
			Latino:     20000.,
		},
	}
	e, err := shp.NewEncoder(TestPopulationShapefile, pop{})
	if err != nil {
		panic(err)
	}
	for _, p := range popData {
		if err = e.Encode(p); err != nil {
			panic(err)
		}
	}
	e.Close()
	f, err := os.Create(strings.TrimRight(TestPopulationShapefile, ".shp") + ".prj")
	if err != nil {
		panic(err)
	}
	if _, err = f.Write([]byte(TestGridSR)); err != nil {
		panic(err)
	}
	f.Close()
}

// WriteTestMortalityShapefile writes out a mortality rate shapefile for testing.
func WriteTestMortalityShapefile() {
	// holder for test mortality data.
	type mortRates struct {
		geom.Polygon
		AllCause, WhNoLMort, BlackMort, AsianMort, NativeMort, LatinoMort float64
	}

	// write out test mortality rate data.
	mortData := []mortRates{
		{
			Polygon: []geom.Path{{
				geom.Point{X: -31.815, Y: -39.4932},
				geom.Point{X: -18.765, Y: -39.4932},
				geom.Point{X: -18.765, Y: -23.1112},
				geom.Point{X: -31.815, Y: -23.1112},
				geom.Point{X: -31.815, Y: -39.4932},
			}},
			AllCause:   800.,
			WhNoLMort:  700.,
			BlackMort:  600.,
			AsianMort:  500.,
			NativeMort: 300.,
			LatinoMort: 400.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -31.815, Y: -23.1112},
				geom.Point{X: -18.765, Y: -23.1112},
				geom.Point{X: -18.765, Y: -6.72923},
				geom.Point{X: -31.815, Y: -6.72923},
				geom.Point{X: -31.815, Y: -23.1112},
			}},
			AllCause:   800.,
			WhNoLMort:  700.,
			BlackMort:  600.,
			AsianMort:  500.,
			NativeMort: 300.,
			LatinoMort: 400.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -31.815, Y: -6.72923},
				geom.Point{X: -18.765, Y: -6.72923},
				geom.Point{X: -18.765, Y: 9.65277},
				geom.Point{X: -31.815, Y: 9.65277},
				geom.Point{X: -31.815, Y: -6.72923},
			}},
			AllCause:   800.,
			WhNoLMort:  700.,
			BlackMort:  600.,
			AsianMort:  500.,
			NativeMort: 300.,
			LatinoMort: 400.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -31.815, Y: 9.65277},
				geom.Point{X: -18.765, Y: 9.65277},
				geom.Point{X: -18.765, Y: 26.03477},
				geom.Point{X: -31.815, Y: 26.03477},
				geom.Point{X: -31.815, Y: 9.65277},
			}},
			AllCause:   800.,
			WhNoLMort:  700.,
			BlackMort:  600.,
			AsianMort:  500.,
			NativeMort: 300.,
			LatinoMort: 400.,
		},
		//Row 2
		{
			Polygon: []geom.Path{{
				geom.Point{X: -18.765, Y: -39.4932},
				geom.Point{X: -5.715, Y: -39.4932},
				geom.Point{X: -5.715, Y: -23.1112},
				geom.Point{X: -18.765, Y: -23.1112},
				geom.Point{X: -18.765, Y: -39.4932},
			}},
			AllCause:   1000.,
			WhNoLMort:  1000.,
			BlackMort:  1000.,
			AsianMort:  1000.,
			NativeMort: 1000.,
			LatinoMort: 1000.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -18.765, Y: -23.1112},
				geom.Point{X: -5.715, Y: -23.1112},
				geom.Point{X: -5.715, Y: -6.72923},
				geom.Point{X: -18.765, Y: -6.72923},
				geom.Point{X: -18.765, Y: -23.1112},
			}},
			AllCause:   1000.,
			WhNoLMort:  1000.,
			BlackMort:  1000.,
			AsianMort:  1000.,
			NativeMort: 1000.,
			LatinoMort: 1000.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -18.765, Y: -6.72923},
				geom.Point{X: -5.715, Y: -6.72923},
				geom.Point{X: -5.715, Y: 9.65277},
				geom.Point{X: -18.765, Y: 9.65277},
				geom.Point{X: -18.765, Y: -6.72923},
			}},
			AllCause:   1000.,
			WhNoLMort:  1000.,
			BlackMort:  1000.,
			AsianMort:  1000.,
			NativeMort: 1000.,
			LatinoMort: 1000.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -18.765, Y: 9.65277},
				geom.Point{X: -5.715, Y: 9.65277},
				geom.Point{X: -5.715, Y: 26.03477},
				geom.Point{X: -18.765, Y: 26.03477},
				geom.Point{X: -18.765, Y: 9.65277},
			}},
			AllCause:   1000.,
			WhNoLMort:  1000.,
			BlackMort:  1000.,
			AsianMort:  1000.,
			NativeMort: 1000.,
			LatinoMort: 1000.,
		},
		//Row 3
		{
			Polygon: []geom.Path{{
				geom.Point{X: -5.715, Y: -39.4932},
				geom.Point{X: 7.335, Y: -39.4932},
				geom.Point{X: 7.335, Y: -23.1112},
				geom.Point{X: -5.715, Y: -23.1112},
				geom.Point{X: -5.715, Y: -39.4932},
			}},
			AllCause:   0.,
			WhNoLMort:  0.,
			BlackMort:  0.,
			AsianMort:  0.,
			NativeMort: 0.,
			LatinoMort: 0.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -5.715, Y: -23.1112},
				geom.Point{X: 7.335, Y: -23.1112},
				geom.Point{X: 7.335, Y: -6.72923},
				geom.Point{X: -5.715, Y: -6.72923},
				geom.Point{X: -5.715, Y: -23.1112},
			}},
			AllCause:   0.,
			WhNoLMort:  0.,
			BlackMort:  0.,
			AsianMort:  0.,
			NativeMort: 0.,
			LatinoMort: 0.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -5.715, Y: -6.72923},
				geom.Point{X: 7.335, Y: -6.72923},
				geom.Point{X: 7.335, Y: 9.65277},
				geom.Point{X: -5.715, Y: 9.65277},
				geom.Point{X: -5.715, Y: -6.72923},
			}},
			AllCause:   0.,
			WhNoLMort:  0.,
			BlackMort:  0.,
			AsianMort:  0.,
			NativeMort: 0.,
			LatinoMort: 0.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: -5.715, Y: 9.65277},
				geom.Point{X: 7.335, Y: 9.65277},
				geom.Point{X: 7.335, Y: 26.03477},
				geom.Point{X: -5.715, Y: 26.03477},
				geom.Point{X: -5.715, Y: 9.65277},
			}},
			AllCause:   0.,
			WhNoLMort:  0.,
			BlackMort:  0.,
			AsianMort:  0.,
			NativeMort: 0.,
			LatinoMort: 0.,
		},
		//Row 4
		{
			Polygon: []geom.Path{{
				geom.Point{X: 7.335, Y: -39.4932},
				geom.Point{X: 20.385, Y: -39.4932},
				geom.Point{X: 20.385, Y: -23.1112},
				geom.Point{X: 7.335, Y: -23.1112},
				geom.Point{X: 7.335, Y: -39.4932},
			}},
			AllCause:   0.,
			WhNoLMort:  0.,
			BlackMort:  0.,
			AsianMort:  0.,
			NativeMort: 0.,
			LatinoMort: 0.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: 7.335, Y: -23.1112},
				geom.Point{X: 20.385, Y: -23.1112},
				geom.Point{X: 20.385, Y: -6.72923},
				geom.Point{X: 7.335, Y: -6.72923},
				geom.Point{X: 7.335, Y: -23.1112},
			}},
			AllCause:   0.,
			WhNoLMort:  0.,
			BlackMort:  0.,
			AsianMort:  0.,
			NativeMort: 0.,
			LatinoMort: 0.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: 7.335, Y: -6.72923},
				geom.Point{X: 20.385, Y: -6.72923},
				geom.Point{X: 20.385, Y: 9.65277},
				geom.Point{X: 7.335, Y: 9.65277},
				geom.Point{X: 7.335, Y: -6.72923},
			}},
			AllCause:   0.,
			WhNoLMort:  0.,
			BlackMort:  0.,
			AsianMort:  0.,
			NativeMort: 0.,
			LatinoMort: 0.,
		},
		{
			Polygon: []geom.Path{{
				geom.Point{X: 7.335, Y: 9.65277},
				geom.Point{X: 20.385, Y: 9.65277},
				geom.Point{X: 20.385, Y: 26.03477},
				geom.Point{X: 7.335, Y: 26.03477},
				geom.Point{X: 7.335, Y: 9.65277},
			}},
			AllCause:   0.,
			WhNoLMort:  0.,
			BlackMort:  0.,
			AsianMort:  0.,
			NativeMort: 0.,
			LatinoMort: 0.,
		},
	}
	e, err := shp.NewEncoder(TestMortalityShapefile, mortRates{})
	if err != nil {
		panic(err)
	}
	for _, m := range mortData {
		if err = e.Encode(m); err != nil {
			panic(err)
		}
	}
	e.Close()
	f, err := os.Create(strings.TrimRight(TestMortalityShapefile, ".shp") + ".prj")
	if err != nil {
		panic(err)
	}
	if _, err = f.Write([]byte(TestGridSR)); err != nil {
		panic(err)
	}
	f.Close()
}

// CreateTestCTMData creates example CTMData for testing. -39.44722, -31.77, 0.09, 0.09
func CreateTestCTMData(load_ctmdata bool) (VarGridConfig, *CTMData) {
	cfg := VarGridConfig{
		VariableGridXo: -39.44722, //-4000,
		VariableGridYo: -31.77,    //-4000,
		VariableGridDx: 13.05,     //4000,
		VariableGridDy: 16.382,    //4000,
		Xnests:         []int{2, 2, 2},
		Ynests:         []int{2, 2, 2},
		HiResLayers:    1,
		//GridProj:            "+proj=lcc +lat_1=33.000000 +lat_2=45.000000 +lat_0=40.000000 +lon_0=-97.000000 +x_0=0 +y_0=0 +a=6370997.000000 +b=6370997.000000 +to_meter=1", "+proj=ob_tran +o_proj=longlat +o_lon_p=87.597031302933 +o_lat_p=31.7583124544932 +lon_0=180 +a=6370997 +units=degrees +no_defs",
		GridProj:            "+proj=longlat",
		PopDensityThreshold: 0.001,
		PopThreshold:        25000,
		PopConcThreshold:    1.0e-7,
		CensusFile:          TestPopulationShapefile,
		CensusPopColumns:    []string{"TotalPop", "WhiteNoLat", "Black", "Native", "Asian", "Latino"},
		PopGridColumn:       "TotalPop",
		MortalityRateFile:   TestMortalityShapefile,
		MortalityRateColumns: map[string]string{
			"AllCause":   "TotalPop",
			"WhNoLMort":  "WhiteNoLat",
			"BlackMort":  "Black",
			"NativeMort": "Native",
			"AsianMort":  "Asian",
			"LatinoMort": "Latino",
		},
	}
	//TR&SB - added switch to load CTM data rather than generating it.
	var cd *CTMData
	var nlayers int
	if !(load_ctmdata) {
		ctmdata := map[string]struct {
			Dims        []string           // netcdf dimensions for this variable
			Description string             // variable description
			Units       string             // variable units
			Data        *sparse.DenseArray // variable data
		}{
			"WindSpeedMinusOnePointFour": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "RMS wind speed^(-1.4)",
				Units:       "(m s-1)^(-1.4)",
			},
			"ParticleDryDep": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Dry deposition velocity for particles",
				Units:       "m s-1",
			},
			"Dz": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Vertical grid size",
				Units:       "m",
			},
			"WAvg": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Annual average z velocity",
				Units:       "m/s",
			},
			"Temperature": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Average Temperature",
				Units:       "K",
			},
			"VOCDryDep": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Dry deposition velocity for VOCs",
				Units:       "m s-1",
			},
			"alt": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Inverse density",
				Units:       "m3 kg-1",
			},
			"UDeviation": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Average deviation from average x velocity",
				Units:       "m/s",
			},
			"gNO": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Average concentration of nitrogen fraction of gaseous NOx",
				Units:       "ug m-3",
			},
			"bVOC": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Average biogenic VOC concentration",
				Units:       "ug m-3",
			},
			"SPartitioning": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Mass fraction of S from SOx in particle {vs. gas} phase",
				Units:       "fraction",
			},
			"Sclass": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Stability parameter",
				Units:       "0=Unstable; 1=Stable",
			},
			"Kzz": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Vertical turbulent diffusivity",
				Units:       "m2 s-1",
			},
			"VDeviation": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Average deviation from average y velocity",
				Units:       "m/s",
			},
			"VAvg": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Annual average y velocity",
				Units:       "m/s",
			},
			"pNH": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Average concentration of nitrogen fraction of particulate ammonium",
				Units:       "ug m-3",
			},
			"SO2DryDep": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Dry deposition velocity for SO2",
				Units:       "m s-1",
			},
			"NOxDryDep": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Dry deposition velocity for NOx",
				Units:       "m s-1",
			},
			"NHPartitioning": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Mass fraction of N from NH3 in particle {vs. gas} phase",
				Units:       "fraction",
			},
			"Kxxyy": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Horizontal eddy diffusion coefficient",
				Units:       "m2 s-1",
			},
			"WindSpeed": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "RMS wind speed",
				Units:       "m s-1",
			},
			"SO2WetDep": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Wet deposition rate constant for SO2 gas",
				Units:       "s-1",
			},
			"LayerHeights": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Height at edge of layer",
				Units:       "m",
			},
			"WindSpeedMinusThird": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "RMS wind speed^(-1/3)",
				Units:       "(m s-1)^(-1/3)",
			},
			"gS": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Average concentration of sulfur fraction of gaseous SOx",
				Units:       "ug m-3",
			},
			"gNH": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Average concentration of nitrogen fraction of gaseous ammonia",
				Units:       "ug m-3",
			},
			"aOrgPartitioning": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Mass fraction of anthropogenic organic matter in particle {vs. gas} phase",
				Units:       "fraction",
			},
			"TotalPM25": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Total PM2.5 concentration",
				Units:       "ug m-3",
			},
			"aVOC": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Average anthropogenic VOC concentration",
				Units:       "ug m-3",
			},
			"NH3DryDep": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Dry deposition velocity for NH3",
				Units:       "m s-1",
			},
			"pS": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Average concentration of sulfur fraction of particulate sulfate",
				Units:       "ug m-3",
			},
			"bOrgPartitioning": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Mass fraction of biogenic organic matter in particle {vs. gas} phase",
				Units:       "fraction",
			},
			"Pblh": {
				Dims:        []string{"y", "x"},
				Data:        sparse.ZerosDense([]int{2, 2}...),
				Description: "Planetary boundary layer height",
				Units:       "m",
			},
			"bSOA": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Average biogenic secondary organic aerosol concentration",
				Units:       "ug m-3",
			},
			"pNO": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Average concentration of nitrogen fraction of particulate NO3",
				Units:       "ug m-3",
			},
			"M2u": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "ACM2 nonlocal upward mixing {Pleim 2007}",
				Units:       "s-1",
			},
			"SO2oxidation": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Rate of SO2 oxidation to SO4 by hydroxyl radical and H2O2",
				Units:       "s-1",
			},
			"OtherGasWetDep": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Wet deposition rate constant for other gases",
				Units:       "s-1",
			},
			"WindSpeedInverse": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "RMS wind speed^(-1)",
				Units:       "(m s-1)^(-1)",
			},
			"NO_NO2partitioning": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Mass fraction of N in NOx that exists as NO.",
				Units:       "fraction",
			},
			"ParticleWetDep": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Wet deposition rate constant for fine particles",
				Units:       "s-1",
			},
			"aSOA": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Average anthropogenic secondary organic aerosol concentration",
				Units:       "ug m-3",
			},
			"NOPartitioning": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Mass fraction of N from NOx in particle {vs. gas} phase",
				Units:       "fraction",
			},
			"M2d": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "ACM2 nonlocal downward mixing {Pleim 2007}",
				Units:       "s-1",
			},
			"S1": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Stability parameter",
				Units:       "?",
			},
			"UAvg": {
				Dims:        []string{"z", "y", "x"},
				Data:        sparse.ZerosDense([]int{10, 2, 2}...),
				Description: "Annual average x velocity",
				Units:       "m/s",
			},
		}

		// Data extracted from a file created by the WRF-Chem preprocessor
		// TODO: This data is missing the staggered information, and should be fixed.
		ctmdata["Dz"].Data.Elements = []float64{55.616737365722656, 55.60960006713867, 55.558589935302734, 55.63689422607422, 80.05044555664062, 80.0164566040039, 79.99669647216797, 80.07341003417969, 105.04544067382812, 104.97071075439453, 104.981201171875, 105.04759216308594, 130.77992248535156, 130.65481567382812, 130.7270050048828, 130.7596435546875, 165.87161254882812, 165.6465301513672, 165.8095703125, 165.7937774658203, 211.19100952148438, 210.81178283691406, 210.9864959716797, 210.94813537597656, 250.37254333496094, 249.94908142089844, 250.0945587158203, 250.04212951660156, 450.19708251953125, 449.4986877441406, 449.73907470703125, 449.6775207519531, 471.3617858886719, 470.7438659667969, 470.8647155761719, 470.8252868652344, 496.6903991699219, 495.80023193359375, 495.7990417480469, 495.63519287109375}
		ctmdata["WAvg"].Data.Elements = []float64{-0.0031770633067935705, -0.03131113201379776, -0.0016643835697323084, -0.022064168006181717, 0.006579460576176643, -0.008502744138240814, 0.0039610546082258224, -0.0008686335058882833, 0.013467869721353054, 0.0042756046168506145, 0.008733110502362251, 0.006209260318428278, 0.01728512905538082, 0.014318481087684631, 0.01065493281930685, 0.008901089429855347, 0.016626669093966484, 0.01868339814245701, 0.00920271035283804, 0.007909671403467655, 0.012121721170842648, 0.017389992251992226, 0.0059646242298185825, 0.005230602342635393, 0.005212655756622553, 0.012798762880265713, 0.0025486990343779325, 0.0049898698925971985, 0.0002958668628707528, 0.010805307887494564, 0.0007378551526926458, 0.010255983099341393, 0.0048363772220909595, 0.008727352134883404, 0.0052039166912436485, 0.016550913453102112, 0.018345091491937637, 0.017693055793642998, 0.0141964852809906, 0.019544612616300583}
		ctmdata["Temperature"].Data.Elements = []float64{281.7968444824219, 281.8467102050781, 281.4716796875, 281.91217041015625, 281.6050109863281, 281.5710144042969, 281.3937683105469, 281.6999816894531, 281.0885314941406, 280.9794006347656, 280.8974609375, 281.1159362792969, 280.289794921875, 280.1217956542969, 280.1652526855469, 280.2762756347656, 279.3124084472656, 279.05828857421875, 279.225341796875, 279.2402648925781, 278.082763671875, 277.7137756347656, 277.8447570800781, 277.838134765625, 276.3756103515625, 276.03662109375, 276.10333251953125, 276.07647705078125, 273.9061584472656, 273.6038818359375, 273.66888427734375, 273.643310546875, 271.0890808105469, 270.8735046386719, 270.8675231933594, 270.8497619628906, 269.2445983886719, 268.8931579589844, 268.8183288574219, 268.7208251953125}
		ctmdata["WindSpeedMinusOnePointFour"].Data.Elements = []float64{0.7428295612335205, 0.5058742165565491, 0.934855043888092, 0.5119651556015015, 0.6354432106018066, 0.5275496244430542, 1.0058796405792236, 0.5661730170249939, 0.636940598487854, 0.6686421036720276, 0.8879344463348389, 0.9787752032279968, 0.8092939257621765, 0.9098353981971741, 1.020829439163208, 0.8001529574394226, 0.7081223726272583, 0.7041523456573486, 0.8658035397529602, 0.739102840423584, 0.5948684215545654, 0.559593677520752, 0.6883410215377808, 0.6498873233795166, 0.5875554084777832, 0.5690711736679077, 1.1264561414718628, 0.8045519590377808, 0.41173169016838074, 0.367610901594162, 0.4526044428348541, 0.2955501675605774, 0.134508416056633, 0.13035288453102112, 0.10757080465555191, 0.11621700972318649, 0.06143762543797493, 0.05990791693329811, 0.061796411871910095, 0.061518121510744095}
		ctmdata["ParticleDryDep"].Data.Elements = []float64{0.0007427706732414663, 0.0007441723137162626, 0.0007016236777417362, 0.0007579709053970873, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
		ctmdata["UDeviation"].Data.Elements = []float64{1.2811052799224854, 1.750707983970642, 1.2451300621032715, 1.491990566253662, 1.7152024507522583, 1.9533977508544922, 1.4286575317382812, 1.2071692943572998, 1.8392454385757446, 1.6557434797286987, 1.5329666137695312, 0.9727783203125, 1.690688133239746, 1.3520514965057373, 1.7273163795471191, 1.1815413236618042, 1.597508430480957, 1.4398936033248901, 1.8474749326705933, 1.610512137413025, 1.444093108177185, 1.4238202571868896, 1.8994537591934204, 1.684370517730713, 1.4248212575912476, 1.3437374830245972, 2.175584554672241, 1.9596269130706787, 2.204925537109375, 1.8786612749099731, 2.3146636486053467, 2.1970036029815674, 2.006307601928711, 1.961408019065857, 1.7362703084945679, 1.693990707397461, 1.9431530237197876, 1.8722081184387207, 1.6691739559173584, 1.6544681787490845}
		ctmdata["gNO"].Data.Elements = []float64{6.137321472167969, 4.699996471405029, 18.284534454345703, 8.36161994934082, 5.439126968383789, 3.9819915294647217, 13.961697578430176, 6.630396842956543, 5.010383605957031, 3.53118896484375, 10.873613357543945, 5.593719005584717, 4.645527362823486, 3.262540102005005, 8.41073226928711, 4.846043109893799, 3.8698558807373047, 2.8284096717834473, 5.868012428283691, 3.5862154960632324, 2.66652512550354, 2.035365104675293, 3.4560422897338867, 2.095470905303955, 1.8920834064483643, 1.1830661296844482, 1.4191572666168213, 0.8098746538162231, 0.934967577457428, 0.6002140045166016, 0.5770926475524902, 0.4145469665527344, 0.3197759687900543, 0.29826974868774414, 0.2326596975326538, 0.23800556361675262, 0.13441583514213562, 0.132355198264122, 0.09638924896717072, 0.11154797673225403}
		ctmdata["bVOC"].Data.Elements = []float64{0.03590616583824158, 0.053436391055583954, 0.051482152193784714, 0.05991194397211075, 0.031461115926504135, 0.04459478706121445, 0.045079413801431656, 0.04948962852358818, 0.027608249336481094, 0.03701688349246979, 0.03900148719549179, 0.0424819141626358, 0.02413773722946644, 0.030795326456427574, 0.03209312632679939, 0.035902272909879684, 0.018761392682790756, 0.022274844348430634, 0.02411295659840107, 0.02881707064807415, 0.013118304312229156, 0.013295278884470463, 0.015942882746458054, 0.01890280283987522, 0.010445994324982166, 0.008873539976775646, 0.00827344786375761, 0.00805355142802, 0.0064531429670751095, 0.00594680430367589, 0.00473548611626029, 0.0048340726643800735, 0.0044800275936722755, 0.0041855378076434135, 0.0037518178578466177, 0.0034344540908932686, 0.002398792887106538, 0.002276006853207946, 0.0019890377297997475, 0.0019814185798168182}
		ctmdata["SPartitioning"].Data.Elements = []float64{0.35624757409095764, 0.4765892028808594, 0.1699211150407791, 0.2879854738712311, 0.3557271957397461, 0.5236513614654541, 0.17515133321285248, 0.3525265157222748, 0.3817720115184784, 0.5159932971000671, 0.21832454204559326, 0.3477080166339874, 0.44939202070236206, 0.5117040872573853, 0.3439377546310425, 0.33935996890068054, 0.4687609374523163, 0.4193888306617737, 0.3161131739616394, 0.37612083554267883, 0.42450302839279175, 0.44209426641464233, 0.4657171368598938, 0.454298198223114, 0.5342852473258972, 0.6493272185325623, 0.5906006097793579, 0.749495804309845, 0.7275176644325256, 0.7907597422599792, 0.804292619228363, 0.8494125008583069, 0.8254836797714233, 0.7912582755088806, 0.8296924829483032, 0.8125612139701843, 0.8629156947135925, 0.7844576239585876, 0.8411360383033752, 0.8015005588531494}
		ctmdata["VOCDryDep"].Data.Elements = []float64{0.0042526074685156345, 0.004445951897650957, 0.0042899735271930695, 0.003794486401602626, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
		ctmdata["alt"].Data.Elements = []float64{0.8143327236175537, 0.8197183012962341, 0.8161125779151917, 0.817907452583313, 0.8204646706581116, 0.8256438374519348, 0.8225621581077576, 0.8240010738372803, 0.8281850814819336, 0.8331775665283203, 0.8303545117378235, 0.8315349817276001, 0.8377386331558228, 0.8425818681716919, 0.8401082158088684, 0.8409796357154846, 0.8500133156776428, 0.8545834422111511, 0.8524419069290161, 0.8530330061912537, 0.8657819628715515, 0.8700550198554993, 0.8677392601966858, 0.8682661652565002, 0.8848150968551636, 0.8892749547958374, 0.8866899013519287, 0.8872032761573792, 0.9148464798927307, 0.9195871353149414, 0.9168713092803955, 0.9174686074256897, 0.9578233361244202, 0.9630182981491089, 0.9599078297615051, 0.960583508014679, 1.0092519521713257, 1.014236569404602, 1.0106980800628662, 1.0111606121063232}
		ctmdata["VAvg"].Data.Elements = []float64{0.3878466784954071, 0.042083218693733215, -0.3559805154800415, -0.5641781091690063, 0.22034627199172974, -0.5800857543945312, -0.3915371298789978, -1.0190788507461548, -0.18792983889579773, -1.0643645524978638, -0.5425596237182617, -1.2878860235214233, -0.5102839469909668, -1.339412808418274, -0.6198769807815552, -1.357798457145691, -0.7527636885643005, -1.4363676309585571, -0.5296856164932251, -1.361451506614685, -0.7092278599739075, -1.2328779697418213, -0.3427528440952301, -1.1049786806106567, -0.8593996167182922, -1.0891087055206299, -0.10885697603225708, -0.6699331402778625, -0.18612952530384064, -0.3922097384929657, 0.42369094491004944, -0.1022019013762474, 1.21445894241333, 0.9678690433502197, 1.261110544204712, 0.7657799124717712, 2.078648567199707, 1.8895270824432373, 2.194763660430908, 1.912973165512085}
		ctmdata["pNH"].Data.Elements = []float64{0.2842402160167694, 0.26102307438850403, 0.3397374153137207, 0.2595134675502777, 0.2826288640499115, 0.2575586438179016, 0.3238135278224945, 0.2557435929775238, 0.2833107113838196, 0.2554378807544708, 0.32004499435424805, 0.2572919726371765, 0.2846950590610504, 0.2558947205543518, 0.3192558288574219, 0.26366671919822693, 0.28077182173728943, 0.25073519349098206, 0.31189262866973877, 0.25861069560050964, 0.244358628988266, 0.21269740164279938, 0.26755237579345703, 0.21742093563079834, 0.19980283081531525, 0.16102035343647003, 0.1629004180431366, 0.13219214975833893, 0.13671869039535522, 0.10027796030044556, 0.0905759260058403, 0.07659262418746948, 0.04779661446809769, 0.045397792011499405, 0.03417603671550751, 0.036950644105672836, 0.018421389162540436, 0.018827475607395172, 0.011643920093774796, 0.014858880080282688}
		ctmdata["SO2DryDep"].Data.Elements = []float64{0.0019483735086396337, 0.0020372653380036354, 0.0019306899048388004, 0.0007309007924050093, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
		ctmdata["NOxDryDep"].Data.Elements = []float64{0.0005772297736257315, 0.0006088464288040996, 0.0005630956147797406, 0.0005128949997015297, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
		ctmdata["NHPartitioning"].Data.Elements = []float64{0.20804473757743835, 0.12692707777023315, 0.24785539507865906, 0.1591263860464096, 0.2646610736846924, 0.19742926955223083, 0.286176472902298, 0.12133932113647461, 0.29103556275367737, 0.24020065367221832, 0.29969221353530884, 0.1824905425310135, 0.3080105483531952, 0.3114902079105377, 0.31419837474823, 0.25988754630088806, 0.35066214203834534, 0.331311970949173, 0.35708558559417725, 0.4001108705997467, 0.5455805659294128, 0.6248354315757751, 0.6231573820114136, 0.6531552672386169, 0.7216328978538513, 0.8237294554710388, 0.8439907431602478, 0.8740569353103638, 0.9286889433860779, 0.9777018427848816, 0.9791176915168762, 0.9791666865348816, 0.9791666865348816, 0.9791666865348816, 0.9791666865348816, 0.9791666865348816, 0.9791666865348816, 0.9791666865348816, 0.9791666865348816, 0.9791666865348816}
		ctmdata["Sclass"].Data.Elements = []float64{0, 0, 0, 0, 0.625, 0.625, 0.6458333134651184, 0.6458333134651184, 0.5416666865348816, 0.4166666567325592, 0.5625, 0.5, 0.2083333283662796, 0.125, 0.4375, 0.1875, 0.2916666567325592, 0.125, 0.3541666567325592, 0.1458333283662796, 0.3333333432674408, 0.1041666641831398, 0.0416666679084301, 0.0625, 0.02083333395421505, 0.0416666679084301, 0, 0, 0.1458333283662796, 0.2083333283662796, 0.1041666641831398, 0.1875, 0.0833333358168602, 0.1875, 0.125, 0.2083333283662796, 0.7708333134651184, 0.7291666865348816, 0.6666666865348816, 0.5625}
		ctmdata["Kzz"].Data.Elements = []float64{0.9650678038597107, 0.7728250026702881, 1.125778317451477, 1.005699634552002, 2.1285557746887207, 1.783519983291626, 2.3763740062713623, 2.1986358165740967, 2.351567029953003, 2.1163718700408936, 2.525076150894165, 2.3255789279937744, 2.364271640777588, 2.227468967437744, 2.4480409622192383, 2.305467367172241, 2.5211145877838135, 2.4369630813598633, 2.4248063564300537, 2.368791103363037, 2.751412868499756, 2.7220346927642822, 2.626596212387085, 2.541729211807251, 2.9064910411834717, 2.9068634510040283, 2.8753504753112793, 2.8457624912261963, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}
		ctmdata["VDeviation"].Data.Elements = []float64{1.6419901847839355, 1.552883505821228, 0.9576117396354675, 1.1282074451446533, 1.4345885515213013, 1.2338838577270508, 1.1126853227615356, 1.1890350580215454, 1.1986955404281616, 1.1466532945632935, 1.1268880367279053, 1.1795437335968018, 0.958810031414032, 1.0351241827011108, 0.9515238404273987, 1.0391275882720947, 0.6394379138946533, 0.883553147315979, 0.7932561039924622, 0.7102439999580383, 0.5965691804885864, 0.6937295794487, 0.6717283725738525, 0.5079746246337891, 0.7716346979141235, 0.7845242023468018, 0.9944789409637451, 0.9164685606956482, 0.9399864077568054, 0.8586839437484741, 0.7474461793899536, 0.8517684936523438, 0.8913931846618652, 1.269628882408142, 1.0596575736999512, 1.0492056608200073, 1.524090051651001, 1.7090855836868286, 1.5799076557159424, 1.610846757888794}
		ctmdata["Kxxyy"].Data.Elements = []float64{0.7575604915618896, 0.9150991439819336, 0.6208891272544861, 0.7524882555007935, 2.3101816177368164, 1.8487164974212646, 2.6796627044677734, 2.2049927711486816, 2.469972610473633, 2.1148242950439453, 2.815176010131836, 2.5295193195343018, 2.37355637550354, 2.24664306640625, 2.6517460346221924, 2.2723093032836914, 2.303529977798462, 2.2230441570281982, 2.2358744144439697, 2.3116273880004883, 2.6886508464813232, 2.6277735233306885, 2.5014443397521973, 2.3839523792266846, 2.8126206398010254, 2.812807083129883, 2.7501752376556396, 2.688507318496704, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3}
		ctmdata["WindSpeed"].Data.Elements = []float64{2.163347005844116, 2.7272017002105713, 1.8843491077423096, 2.5613532066345215, 2.466365337371826, 2.60990047454834, 1.908056616783142, 2.216294288635254, 2.3336946964263916, 2.297699451446533, 1.7950994968414307, 1.906171202659607, 2.100137948989868, 2.093118190765381, 1.8263593912124634, 1.8500407934188843, 2.0755155086517334, 2.1254751682281494, 1.959427833557129, 1.9774097204208374, 1.9850538969039917, 1.9250260591506958, 2.095578193664551, 2.096379041671753, 1.9812132120132446, 2.021040916442871, 2.520017623901367, 2.5633764266967773, 3.3489553928375244, 3.2535643577575684, 3.8730413913726807, 3.8324387073516846, 5.816560745239258, 5.697556972503662, 5.885941028594971, 5.655944347381592, 7.861310005187988, 7.990108966827393, 7.697640419006348, 7.758955955505371}
		ctmdata["gS"].Data.Elements = []float64{0.2634888291358948, 0.21905501186847687, 0.9337737560272217, 0.41840294003486633, 0.25802916288375854, 0.2025141716003418, 0.7405915856361389, 0.35498514771461487, 0.25728145241737366, 0.19435034692287445, 0.6016327142715454, 0.31845349073410034, 0.25473979115486145, 0.19254085421562195, 0.4821411371231079, 0.29016828536987305, 0.22725290060043335, 0.17923219501972198, 0.3476753532886505, 0.22313575446605682, 0.16836073994636536, 0.13953252136707306, 0.21109969913959503, 0.13692250847816467, 0.1287701427936554, 0.09053260833024979, 0.09876422584056854, 0.06713377684354782, 0.0660431832075119, 0.05449407547712326, 0.04652734845876694, 0.04346112906932831, 0.0412360355257988, 0.04184114560484886, 0.03807821124792099, 0.03864813223481178, 0.043244291096925735, 0.04363596439361572, 0.04144774749875069, 0.0420348159968853}
		ctmdata["gNH"].Data.Elements = []float64{1.1363224983215332, 0.8005860447883606, 2.091249942779541, 0.9906166195869446, 7.3356591868400574, 0.5696154832839966, 1.5564550161361694, 0.7306267023086548, 0.6650115251541138, 0.42707446217536926, 1.1549714803695679, 0.5827423334121704, 0.5340535044670105, 0.3396942615509033, 0.8416059613227844, 0.47406280040740967, 0.3696293234825134, 0.23472052812576294, 0.5016713738441467, 0.30318233370780945, 0.1975899189710617, 0.1367780864238739, 0.1922084391117096, 0.11465863883495331, 0.1045338436961174, 0.039635926485061646, 0.043673817068338394, 0.018766969442367554, 0.0036989825312048197, 0.0001889436534838751, 6.824755018897122e-06, 7.628781233028764e-14, 7.307202802911417e-14, 7.267786632930706e-14, 7.291474417520441e-14, 7.28625195118085e-14, 6.935424456955275e-14, 6.901303936960798e-14, 6.925577190723675e-14, 6.922344235370595e-14}
		ctmdata["aOrgPartitioning"].Data.Elements = []float64{0.00594042195007205, 0.005416125059127808, 0.0030360654927790165, 0.0018056770786643028, 0.0058972700498998165, 0.004511366598308086, 0.010038474574685097, 0.0014150284696370363, 0.007893022149801254, 0.003134393598884344, 0.01587115414440632, 0.001612423686310649, 0.0270244088023901, 0.006842840928584337, 0.02319810353219509, 0.009840000420808792, 0.029563654214143753, 0.02937074564397335, 0.0035459278151392937, 0.0065204850398004055, 0.005085194483399391, 0.005293807480484247, 0.004529155790805817, 0.011838467791676521, 0.016023138538002968, 0.030372390523552895, 0.00795662496238947, 0.01085809338837862, 0.028785619884729385, 0.007128694094717503, 0.010830596089363098, 0.020842688158154488, 0.007270642556250095, 0.007605940569192171, 0.007371330633759499, 0.011374425143003464, 0.022471880540251732, 0.02131063863635063, 0.03438222035765648, 0.05923999845981598}
		ctmdata["TotalPM25"].Data.Elements = []float64{4.907700538635254, 4.25741720199585, 10.347429275512695, 5.3623223304748535, 4.673975944519043, 3.944087505340576, 8.390427589416504, 4.814750671386719, 4.526025772094727, 3.746352434158325, 7.087573051452637, 4.509608745574951, 4.385721206665039, 3.634274482727051, 6.099825382232666, 4.28328275680542, 4.027258396148682, 3.39872670173645, 5.009143352508545, 3.7652153968811035, 3.3942463397979736, 2.9605259895324707, 3.7908523082733154, 2.9806408882141113, 2.9653866291046143, 2.4340357780456543, 2.572143793106079, 2.189840316772461, 2.291722297668457, 1.9475696086883545, 1.9517537355422974, 1.7944254875183105, 1.5315279960632324, 1.4643031358718872, 1.4228588342666626, 1.4233944416046143, 1.0344351530075073, 1.0030040740966797, 0.9847849607467651, 0.9959700703620911}
		ctmdata["aVOC"].Data.Elements = []float64{38.87417221069336, 37.30889892578125, 90.528076171875, 69.00169372558594, 35.89360046386719, 32.94846725463867, 72.40897369384766, 53.24059295654297, 33.76029968261719, 29.89646339416504, 58.896583557128906, 44.11149215698242, 31.78070640563965, 27.97871208190918, 47.88663864135742, 37.819522857666016, 27.568506240844727, 24.56879234313965, 35.51859664916992, 29.220561981201172, 20.18500328063965, 18.46787452697754, 22.517818450927734, 18.731372833251953, 15.217564582824707, 12.044897079467773, 12.207158088684082, 9.69643783569336, 10.008086204528809, 7.978878498077393, 7.933547496795654, 6.906917095184326, 5.5251054763793945, 5.28988790512085, 4.948215484619141, 4.9634857177734375, 3.6040360927581787, 3.56182861328125, 3.346196413040161, 3.473365068435669}
		ctmdata["SO2WetDep"].Data.Elements = []float64{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.2456439559421484e-16, 0, 2.628669343014787e-15, 9.334889649058363e-16, 1.2176135886230247e-14, 1.849259440511123e-14, 2.9082242663331856e-14, 2.2905488676573313e-14, 2.5171894733541644e-14, 4.641805264270736e-14, 8.353235637914569e-15, 2.5565649857534067e-14, 9.540490379861873e-15, 2.0349985192509405e-14, 0, 0, 0, 0}
		ctmdata["LayerHeights"].Data.Elements = []float64{0, 0, 0, 0, 55.616737365722656, 55.60960006713867, 55.558589935302734, 55.63689422607422, 135.66717529296875, 135.6260528564453, 135.55528259277344, 135.71029663085938, 240.71261596679688, 240.59677124023438, 240.53648376464844, 240.75790405273438, 371.4925537109375, 371.2515869140625, 371.26348876953125, 371.5175476074219, 537.3641357421875, 536.8981323242188, 537.0730590820312, 537.3113403320312, 748.55517578125, 747.7098999023438, 748.0595703125, 748.2594604492188, 998.927734375, 997.6589965820312, 998.1541137695312, 998.3015747070312, 1449.124755859375, 1447.15771484375, 1447.8931884765625, 1447.9791259765625, 1920.486572265625, 1917.9014892578125, 1918.7579345703125, 1918.804443359375}
		ctmdata["WindSpeedMinusThird"].Data.Elements = []float64{0.8542997241020203, 0.7842972874641418, 0.9090532660484314, 0.7905508279800415, 0.8236940503120422, 0.7991727590560913, 0.9096444249153137, 0.8206126093864441, 0.8376582264900208, 0.8331751823425293, 0.9196218848228455, 0.8824751973152161, 0.8804322481155396, 0.8712295889854431, 0.9322598576545715, 0.8896374106407166, 0.8621595501899719, 0.8406710028648376, 0.9091718196868896, 0.8771911263465881, 0.8525558710098267, 0.8472396731376648, 0.8702022433280945, 0.8548160791397095, 0.8503544330596924, 0.8395950198173523, 0.8561369180679321, 0.8324417471885681, 0.7381977438926697, 0.732713520526886, 0.7255184650421143, 0.693711519241333, 0.5830445289611816, 0.5849269032478333, 0.5706101655960083, 0.5798718929290771, 0.5093699097633362, 0.5065743327140808, 0.5113828182220459, 0.5105403661727905}
		ctmdata["bOrgPartitioning"].Data.Elements = []float64{0.2325323075056076, 0.1578335016965866, 0.2650301158428192, 0.17489692568778992, 0.22279049456119537, 0.15315434336662292, 0.24865323305130005, 0.1898593306541443, 0.20505140721797943, 0.19346097111701965, 0.26499444246292114, 0.22019386291503906, 0.24375058710575104, 0.19083625078201294, 0.3246062099933624, 0.20013076066970825, 0.2295886129140854, 0.185774028301239, 0.21336442232131958, 0.20717763900756836, 0.22818411886692047, 0.20130540430545807, 0.18484775722026825, 0.1869724690914154, 0.2271156907081604, 0.23805075883865356, 0.2365027815103531, 0.2087409496307373, 0.24162143468856812, 0.20563393831253052, 0.1848113089799881, 0.1625615954399109, 0.11139049381017685, 0.13769954442977905, 0.1571400761604309, 0.1736593246459961, 0.10154013335704803, 0.10977180302143097, 0.10533204674720764, 0.1320185363292694}
		ctmdata["Pblh"].Data.Elements = []float64{213.3524627685547, 238.55056762695312, 213.1160430908203, 233.00067138671875}
		ctmdata["NH3DryDep"].Data.Elements = []float64{0.0006914465920999646, 0.0007092549349181354, 0.0006794998189434409, 0.0003023374010808766, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
		ctmdata["pS"].Data.Elements = []float64{0.284243106842041, 0.2770530581474304, 0.30431225895881653, 0.27223989367485046, 0.2839638590812683, 0.2744355797767639, 0.2966228425502777, 0.2696992754936218, 0.2809333801269531, 0.2711707353591919, 0.29102984070777893, 0.2681629955768585, 0.27604663372039795, 0.26832398772239685, 0.2864548861980438, 0.26791390776634216, 0.2736201286315918, 0.26786676049232483, 0.28727754950523376, 0.27330756187438965, 0.2844190001487732, 0.27556902170181274, 0.29721346497535706, 0.2844560742378235, 0.2986434996128082, 0.2896102964878082, 0.2983483076095581, 0.29068583250045776, 0.3208744525909424, 0.30125153064727783, 0.3066104054450989, 0.2917834520339966, 0.28072044253349304, 0.26913735270500183, 0.2669110894203186, 0.265718013048172, 0.21148428320884705, 0.20438693463802338, 0.20336949825286865, 0.2031562775373459}
		ctmdata["M2u"].Data.Elements = []float64{5.7540772104403004e-05, 0.00011617777636274695, 2.2268945031100884e-05, 3.756620208150707e-05, 5.7540772104403004e-05, 0.00011617777636274695, 2.2268945031100884e-05, 3.756620208150707e-05, 2.95566969725769e-05, 3.502612526062876e-05, 2.2268945031100884e-05, 1.1333466318319552e-05, 1.058194538927637e-05, 1.0673002179828472e-05, 1.495724791311659e-05, 1.1333466318319552e-05, 4.433019967109431e-06, 4.947029992763419e-06, 1.495724791311659e-05, 5.1041238293692e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
		ctmdata["SO2oxidation"].Data.Elements = []float64{1.0774373038202611e-07, 1.3951567723324843e-07, 6.831863430534213e-08, 1.2388677816943527e-07, 1.0202472111586758e-07, 1.3523620623345778e-07, 6.534514085387855e-08, 1.1084567574926041e-07, 1.020713895627523e-07, 1.3866753079128102e-07, 6.780697958674864e-08, 1.0870659394868198e-07, 1.0652390614040996e-07, 1.4826964900294115e-07, 7.382872979633248e-08, 1.1307518832381902e-07, 1.327407801454683e-07, 1.8692415437726595e-07, 9.605693662706472e-08, 1.3681913912932941e-07, 2.729773598275642e-07, 2.9107002319506137e-07, 2.0556309721087018e-07, 2.0854740512277203e-07, 4.0155558167498384e-07, 4.475469381759467e-07, 3.852896952594165e-07, 4.407461631217302e-07, 5.712989263884083e-07, 5.519483465832309e-07, 1.7857456668934901e-06, 1.5178520698100328e-06, 1.196515086121508e-06, 1.1856199080284568e-06, 1.1644171991065377e-06, 1.1659409437925206e-06, 6.919872816979478e-07, 6.666737704108527e-07, 6.51064908652188e-07, 6.511202741421585e-07}
		ctmdata["OtherGasWetDep"].Data.Elements = []float64{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7.4854795002101e-16, 0, 8.762231425726939e-15, 3.1116299536055334e-15, 4.0587118491390226e-14, 6.164198247974803e-14, 9.694080661901833e-14, 7.635163118066557e-14, 8.39063112609631e-14, 1.547268444011124e-13, 2.78441176636713e-14, 8.521883285844689e-14, 3.180163516422821e-14, 6.783328510440861e-14, 0, 0, 0, 0}
		ctmdata["bSOA"].Data.Elements = []float64{0.0011581587605178356, 0.000992136774584651, 0.002221859758719802, 0.0012743411352857947, 0.0011546823661774397, 0.0009657172486186028, 0.0020205071195960045, 0.001230109017342329, 0.0011784618254750967, 0.0009719033841975033, 0.0019361464073881507, 0.0012480099685490131, 0.0012439475394785404, 0.0010027826065197587, 0.0019435686990618706, 0.0013141032541170716, 0.0013006493682041764, 0.0009734203922562301, 0.0019110259599983692, 0.0013258332619443536, 0.001143123721703887, 0.0008504534489475191, 0.001599263516254723, 0.0010998825309798121, 0.0008840354275889695, 0.000671887188218534, 0.0009638185147196054, 0.0007680259295739233, 0.0004998059594072402, 0.00042432971531525254, 0.00043949956307187676, 0.0004122521204408258, 0.00021422718418762088, 0.00021404004655778408, 0.00018060504226014018, 0.0001923404197441414, 7.241401908686385e-05, 7.528419519076124e-05, 6.001858855597675e-05, 6.772561027901247e-05}
		ctmdata["pNO"].Data.Elements = []float64{0.03590570390224457, 0.018970012664794922, 0.07387321442365646, 0.02166839875280857, 0.03453795611858368, 0.018142065033316612, 0.06466609984636307, 0.020116569474339485, 0.03786734491586685, 0.019568972289562225, 0.06578326225280762, 0.02300654910504818, 0.04352101683616638, 0.023323234170675278, 0.07103673368692398, 0.02972385659813881, 0.0435439869761467, 0.023871595039963722, 0.06975264102220535, 0.03307438641786575, 0.03700215741991997, 0.02524322271347046, 0.050182048231363297, 0.02205708995461464, 0.04037543758749962, 0.01847105659544468, 0.0196200180798769, 0.010634751990437508, 0.01481037400662899, 0.003672539023682475, 0.0023068245500326157, 0.0016349966172128916, 0.0006205170648172498, 0.00038912950549274683, 0.0005975642125122249, 0.0005193596589379013, 8.650866948300973e-05, 5.208631409914233e-05, 0.00011906155123142526, 0.00011644966434687376}
		ctmdata["aSOA"].Data.Elements = []float64{0.13039030134677887, 0.11563628166913986, 0.22614403069019318, 0.14091438055038452, 0.1277293711900711, 0.11026536673307419, 0.19979557394981384, 0.1331716626882553, 0.1263304352760315, 0.10758806765079498, 0.1793670505285263, 0.12926062941551208, 0.12510785460472107, 0.10693610459566116, 0.16245737671852112, 0.12507545948028564, 0.11649894714355469, 0.10026714205741882, 0.1382245421409607, 0.10982552170753479, 0.09399762004613876, 0.08427876979112625, 0.10569468885660172, 0.08468136936426163, 0.07855755090713501, 0.06759165227413177, 0.07104936987161636, 0.05824333056807518, 0.0625729039311409, 0.05500199645757675, 0.05291100591421127, 0.0478048212826252, 0.039398252964019775, 0.03919265791773796, 0.03528725728392601, 0.03726901113986969, 0.02428617887198925, 0.024222973734140396, 0.023149261251091957, 0.024320900440216064}
		ctmdata["NOPartitioning"].Data.Elements = []float64{0.03755469620227814, 0.008940469473600388, 0.029697928577661514, 0.020699674263596535, 0.018331794068217278, 0.014127560891211033, 0.029086602851748466, 0.0077688442543148994, 0.029589436948299408, 0.009707137942314148, 0.015623725950717926, 0.008958968333899975, 0.056342996656894684, 0.018050367012619972, 0.021349797025322914, 0.044204212725162506, 0.05978367105126381, 0.04683312773704529, 0.03634558245539665, 0.03328275680541992, 0.0314299650490284, 0.056021589785814285, 0.043720345944166183, 0.03292354196310043, 0.09827729314565659, 0.052806224673986435, 0.04945669323205948, 0.09483204036951065, 0.04705972597002983, 0.039231300354003906, 0.031066715717315674, 0.02624012902379036, 0.013843486085534096, 0.004339892882853746, 0.0034604202955961227, 0.0030168064404278994, 0.00025264229043386877, 0.00038729843799956143, 0.0015666232211515307, 0.0055247205309569836}
		ctmdata["M2d"].Data.Elements = []float64{0.00023428934218827635, 0.00038927432615309954, 0.00017624394968152046, 0.00015487683413084596, 0.00012279980001039803, 0.00018979581363964826, 0.00010693733202060685, 8.151019574142992e-05, 4.9730995669960976e-05, 5.611712549580261e-05, 6.451813533203676e-05, 3.349667531438172e-05, 1.6204461644520052e-05, 1.69449358509155e-05, 3.392849612282589e-05, 1.7805126844905317e-05, 4.433019967109431e-06, 4.947029992763419e-06, 1.495724791311659e-05, 5.1041238293692e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
		ctmdata["S1"].Data.Elements = []float64{0, 0, 0, 0, 2.9732120310654864e-05, 2.43374324782053e-05, 3.705010749399662e-05, 2.8398371796356514e-05, 1.688007978373207e-05, 1.3518683772417717e-05, 1.7779218978830613e-05, 1.3869238500774372e-05, 1.1639571312116459e-05, 9.61807290877914e-06, 1.3893534742237534e-05, 1.024402354232734e-05, 1.2529078958323225e-05, 1.0150060916203074e-05, 1.353778407064965e-05, 1.0917246981989592e-05, 1.2855670320277568e-05, 1.0336762898077723e-05, 9.572671842761338e-06, 9.103485353989527e-06, 8.982809958979487e-06, 9.45913780014962e-06, 8.373491255042609e-06, 8.023001100809779e-06, 1.3160735761630349e-05, 1.36547205329407e-05, 1.3641458281199448e-05, 1.3657489034812897e-05, 1.3380425116338301e-05, 1.4064138667890802e-05, 1.3495051462086849e-05, 1.3554791621572804e-05, 2.2211648683878593e-05, 2.1147137886146083e-05, 2.061365739791654e-05, 1.9986231563962065e-05}
		ctmdata["UAvg"].Data.Elements = []float64{-0.041519902646541595, -0.9808820486068726, 0.06470483541488647, -0.7718557119369507, -0.3607710301876068, -0.9662059545516968, 0.1872982680797577, -0.38807952404022217, -0.5591860413551331, -0.6723234057426453, 0.3609139323234558, 0.18215824663639069, -0.6170716881752014, -0.44226333498954773, 0.47974804043769836, 0.5847874283790588, -0.6786242127418518, -0.4665543735027313, 0.4387609660625458, 0.6910924315452576, -0.4648977220058441, -0.4054996967315674, 0.6767528653144836, 0.7725197672843933, 0.8310025334358215, 0.5767632126808167, 1.6348024606704712, 1.6891616582870483, 3.1902570724487305, 2.6577134132385254, 3.724231719970703, 3.490123987197876, 5.821722030639648, 5.228752136230469, 5.773016929626465, 5.378808498382568, 7.232807159423828, 7.238331317901611, 7.092085361480713, 7.030161380767822}
		ctmdata["WindSpeedInverse"].Data.Elements = []float64{0.732271134853363, 0.5603066682815552, 0.8760741353034973, 0.5667445659637451, 0.657837986946106, 0.5875391960144043, 0.9011663198471069, 0.620093822479248, 0.6740943789482117, 0.681207537651062, 0.8681759834289551, 0.8433297276496887, 0.794369637966156, 0.815617024898529, 0.9353957176208496, 0.7946033477783203, 0.7256165146827698, 0.6956785321235657, 0.8492530584335327, 0.7572640776634216, 0.6654402017593384, 0.6422377228736877, 0.7287707924842834, 0.6930571794509888, 0.6594492197036743, 0.6392582058906555, 0.8657473921775818, 0.7384691834449768, 0.47879090905189514, 0.45145541429519653, 0.4920858144760132, 0.38517093658447266, 0.22066578269004822, 0.2190883904695511, 0.19627341628074646, 0.20691488683223724, 0.13475628197193146, 0.13244254887104034, 0.1357078105211258, 0.1351943463087082}
		ctmdata["NO_NO2partitioning"].Data.Elements = []float64{0.8928555846214294, 0.8909785151481628, 0.641761302947998, 0.8497637510299683, 0.9013280272483826, 0.9032465815544128, 0.7803301811218262, 0.8758372068405151, 0.902558445930481, 0.8780083656311035, 0.8402559161186218, 0.8790870904922485, 0.8564733266830444, 0.8600209951400757, 0.873668909072876, 0.8516348600387573, 0.8567376732826233, 0.9022117257118225, 0.8590437769889832, 0.901540994644165, 0.908869743347168, 0.8370416164398193, 0.9239614009857178, 0.9153663516044617, 0.9301338195800781, 0.9184640645980835, 0.9023407101631165, 0.882240891456604, 0.8637320399284363, 0.8614365458488464, 0.9050966501235962, 0.8893399834632874, 0.8403955101966858, 0.8611547350883484, 0.8262611627578735, 0.8487136960029602, 0.8544286489486694, 0.8645439743995667, 0.8191149234771729, 0.8504053950309753}
		ctmdata["ParticleWetDep"].Data.Elements = []float64{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.2850170500475784e-11, 0, 4.5631501355281046e-10, 1.6177180428567084e-10, 2.1113284454088443e-09, 3.2059230825609575e-09, 4.5525401226598206e-09, 4.860858382471633e-09, 3.966416084466573e-09, 7.3075998585636626e-09, 1.3649917862679217e-09, 4.172174605798773e-09, 1.6363523869245e-09, 3.835154860354351e-09, 0, 0, 0, 0}
		cd = &CTMData{
			Data: ctmdata,
			xo:   -12000,
			yo:   -12000,
			dx:   12000,
			dy:   12000,
			nx:   2,
			ny:   2,
		}
		nlayers = 10
	} else {
		loadFileName := "cmd/inmap/testdata/preproc/inmapData_GEMMACH_golden.nc"
		f2, err := os.Open(loadFileName)
		if err != nil {
			print(err)
		}
		cd, err = cfg.LoadCTMData(f2)
		if err != nil {
			print(err)
		}
		nlayers = 84
	}

	cd.makeCTMgrid(nlayers)

	return cfg, cd
}

// VarGridTestData returns some test data for variable grid generation.
func VarGridTestData() (*VarGridConfig, *CTMData, *Population, PopIndices, *MortalityRates, MortIndices) {

	WriteTestPopShapefile()
	WriteTestMortalityShapefile()
	load_ctmdata := true
	cfg, data := CreateTestCTMData(load_ctmdata)

	population, popIndices, mortalityRates, mortIndices, err := cfg.LoadPopMort()
	if err != nil {
		panic(err)
	}

	for _, fname := range []string{TestPopulationShapefile, TestMortalityShapefile} {
		DeleteShapefile(fname)
	}

	return &cfg, data, population, popIndices, mortalityRates, mortIndices
}

// DeleteShapefile deletes the named shapefile.
func DeleteShapefile(fname string) error {
	for _, ext := range []string{".dbf", ".prj", ".shp", ".shx"} {
		if err := os.Remove(strings.TrimSuffix(fname, ".shp") + ext); err != nil {
			return err
		}
	}
	return nil
}
