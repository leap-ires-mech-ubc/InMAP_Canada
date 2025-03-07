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
	"math"
	"os"
	"strconv"
	"strings"

	//"math"
	"time"

	"github.com/ctessum/atmos/seinfeld"
	"github.com/ctessum/atmos/wesely1989"
	"github.com/ctessum/cdf"
	"github.com/ctessum/geom"
	"github.com/ctessum/geom/index/rtree"

	"github.com/ctessum/sparse"
	//TR added 20220422 for aSOA/bSOA partition coefficient estimation. Could also copy to keep internal.
	//"github.com/sajari/regression"
)

//TR - update to GEMMACH variables
// GEM-MACH variables currently used:
/* hc5,hc8,olt,oli,tol,xyl,csl,cvasoa1,cvasoa2,cvasoa3,cvasoa4,iso,api,sesq,lim,
cvbsoa1,cvbsoa2,cvbsoa3,cvbsoa4,asoa1i,asoa1j,asoa2i,asoa2j,asoa3i,asoa3j,asoa4i,
asoa4j,bsoa1i,bsoa1j,bsoa2i,bsoa2j,bsoa3i,bsoa3j,bsoa4i,bsoa4j,no,no2,no3ai,no3aj,
so2,sulf,so4ai,so4aj,nh3,nh4ai,nh4aj,PM2_5_DRY,U,V,W,PBLH,PH,PHB,HFX,UST,T,
PB,P,ho,h2o2,LU_INDEX,QRAIN,CLDFRA,QCLOUD,ALT,SWDOWN,GLW */

//const wrfFormat = "2006-01-02_15_04_05"
//20220907 Changed to GEMMACH time format. Pulled from test netcdf file provided by ECCC
//const gemFormat = "2006-01-02T15:04:05.000000000"
//TR 20221115 - changed to file naming convention. Still having prboems with file reader - need to figure out later perhaps
const (
	gemFormat     = "2006-01-02_15_04_05" //"2019010100_000"
	gemChemFormat = "20060102.150000"
)

var Nx, Ny, Nz int

// GEMMACH is an InMAP preprocessor for GEM-MACH output.
//recordDelta, fileDelta time.Duration
type GEMMACH struct {
	aVOC, bVOC, tSOA, aSOA, bSOA, nox, no, no2, pNO, sox, pS, nh3, pNH, totalPM25, ho, h2o2 map[string]float64

	gemnoChemHour bool

	start, end time.Time

	gem_out, gem_rdps string

	//gem_geophy string

	recordDelta, fileDelta time.Duration

	msgChan chan string

	//SB 20221127
	//start Copy-Adapt-VF-v0: copied variables needed for vegetation fraction file processing from geoschem.go
	landUse *sparse.DenseArray

	nz int

	dx, dy float64

	xCenters, yCenters []float64
	//end Copy-Adapt-VF-v0

}

// NewGEMMACH initializes a GEM-MACH preprocessor from the given
// configuration information.
// gem_out is the location of GEM-MACH output files.
// [DATE] should be used as a wild card for the simulation date.
// startDate and endDate are the dates of the beginning and end of the
// simulation, respectively, in the format "YYYYMMDD".
// If msgChan is not nil, status messages will be sent to it.
//func NewGEMMACH(gem_out, startDate, endDate string, noChemHour bool, msgChan chan string) (*GEMMACH, error) {

//SB 20221127:
//start Copy-Adapt-VF-v0
//func NewGEOSChem(GEOSA1, GEOSA3Cld, GEOSA3Dyn, GEOSI3, GEOSA3MstE, GEOSApBp, GEOSChemOut, OlsonLandMap, startDate, endDate string, dash bool, chemRecordStr, chemFileStr string, noChemHour bool, msgChan chan string) (*GEOSChem, error) {
func NewGEMMACH(gem_out, gem_geophy, gem_rdps, startDate, endDate string, gemnoChemHour bool, msgChan chan string) (*GEMMACH, error) {
	//end Copy-Adapt-VF-v0
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
			"TA2": 0.05, "TA3": 0.05, "TTOL": 1, "TALD": 1, "TARO": 1, "TC38": 1, "TCRE": 1,
			"TDIA": 1, "TETH": 1, "THCH": 1, "TMGL": 1,
		},
		//"TTOL": 1,"TA3": 1,"TA2": 1,

		//TR - Need to adapt to GEMMACH aVOC variables - to be confirmed by ECCC
		bVOC: map[string]float64{
			"TA2": 0.95, "TA3": 0.95, "TH22": 1, "TISO": 1,
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
		aSOA: map[string]float64{"TOC1": 1.},
		// VBS SOA species (biogenic only) [μg/kg dry air].
		//bSOA: map[string]float64{"bSOA": SOAest(tSOA, aVOC, bVOC, false)},
		bSOA: map[string]float64{"TOC1": 1.},
		// NOx is RACM NOx species. We are only interested in the mass
		// of Nitrogen, rather than the mass of the whole molecule, so
		// we use the molecular weight of Nitrogen.
		nox: map[string]float64{"TNO": mwN / mwNO,
			"TNO2": mwN / mwNOx},
		// pNO is the Nitrogen fraction of MADE particulate
		// NO species [μg/kg dry air].
		pNO: map[string]float64{"TNI1": mwN / mwNO3,
			"TNO3": mwN / mwNO3,
			"THN3": mwN / mwHN3},
		// SOx is the RACM SOx species. We are only interested in the mass
		// of Sulfur, rather than the mass of the whole molecule, so
		// we use the molecular weight of Sulfur.
		//TR/SB - converted from PPB to PPM by dividing by 1000
		sox: map[string]float64{"S2": ppmvToUgKg(mwS) / 1000.0},
		// pS is the Sulfur fraction of the MADE particulate
		// Sulfur species [μg/kg dry air].
		pS: map[string]float64{"TSU1": mwS / mwSO4},
		// NH3 is ammonia. We are only interested in the mass
		// of Nitrogen, rather than the mass of the whole molecule, so
		// we use the molecular weight of Nitrogen.
		nh3: map[string]float64{"TNH3": mwN / mwNH3},
		// pNH is the Nitrogen fraction of the MADE particulate
		// ammonia species [μg/kg].
		pNH: map[string]float64{"TAM1": mwN / mwNH4},
		// totalPM25 is total mass of PM2.5  [μg/m3].
		totalPM25: map[string]float64{"AF": 1.},
		// Hydroxy radical is concentratoon of the hydroxy radical.
		//Convert from ug/kg to ppmV for consistency
		ho: map[string]float64{"TOH": UgKgToppmv(mwOH)},
		//hydrogen peroxide is concentratoon of H202.
		//Convert from ug/kg to ppmV for consistency
		h2o2: map[string]float64{"TH22": UgKgToppmv(mwH2O2)},

		gem_out:       gem_out,
		gem_rdps:      gem_rdps,
		msgChan:       msgChan,
		gemnoChemHour: gemnoChemHour,
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
	w.fileDelta, err = time.ParseDuration("1h")
	if err != nil {
		return nil, fmt.Errorf("inmap: GEM-MACH preprocessor fileDelta: %v", err)
	}

	//SB 20221127:
	//start Copy-Adapt-VF-v0
	file, err := os.Open(gem_geophy)
	if err != nil {
		return nil, err
	}
	//commenting below defer file.Close() since will keep uncommented the last defer file.close in the following copy-adapt-VF-v0
	//if does not work, simply uncomment all the commented sections here
	//defer file.Close()
	cfile, err := cdf.Open(file)
	if err != nil {
		return nil, fmt.Errorf("inmap: gem_geophy land use file: %v", err)
	}

	//commenting below lines to try a single os.Open, cdf.Open, and defer file.Close() in the following sections of copy-adapt-VF-v0
	//if does not work, simply uncomment all the commented sections here
	//file, err := os.Open(Gem_geophy)
	//if err != nil {
	//	return nil, err
	//}
	//defer file.Close()
	//cfile, err := cdf.Open(file)
	//if err != nil {
	//	return nil, fmt.Errorf("inmap: Gem_geophy land use file: %v", err)
	//}

	//w.dx, err = w.DX()
	w.dx, err = DX(cfile)
	if err != nil {
		return nil, err
	}
	//print(w.dx)

	//file, err := os.Open(Gem_geophy)
	//if err != nil {
	//	return nil, err
	//}
	//defer file.Close()
	//cfile, err := cdf.Open(file)
	//if err != nil {
	//	return nil, fmt.Errorf("inmap: Gem_geophy land use file: %v", err)
	//}
	//w.dy, err = w.DY()
	w.dy, err = DY(cfile)
	if err != nil {
		return nil, err
	}
	//print(w.dy)

	//file, err := os.Open(Gem_geophy)
	//if err != nil {
	//	return nil, err
	//}
	//defer file.Close()
	//cfile, err := cdf.Open(file)
	//if err != nil {
	//	return nil, fmt.Errorf("inmap: Gem_geophy land use file: %v", err)
	//}
	//w.xCenters, err = w.xCenters()
	w.xCenters, err = xCenters(cfile)
	if err != nil {
		return nil, err
	}
	//print(w.xCenters)
	//file, err := os.Open(Gem_geophy)
	//if err != nil {
	//	return nil, err
	//}
	//defer file.Close()
	//cfile, err := cdf.Open(file)
	//if err != nil {
	//	return nil, fmt.Errorf("inmap: Gem_geophy land use file: %v", err)
	//}
	//w.yCenters, err = w.yCenters()
	w.yCenters, err = yCenters(cfile)
	if err != nil {
		return nil, err
	}
	//print(w.yCenters)
	//file, err := os.Open(Gem_geophy)
	//if err != nil {
	//	return nil, err
	//}
	//defer file.Close()
	//cfile, err := cdf.Open(file)
	//if err != nil {
	//	return nil, fmt.Errorf("inmap: Gem_geophy land use file: %v", err)
	//}
	w.landUse, err = w.largestLandUse(cfile)
	if err != nil {
		return nil, err
	}
	//end Copy-Adapt-VF-v0

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
//TR: Syntax for these methods since I find them a bit confusing. Basically, a method is a function for a struct (or similar)
//These functions set up w as a "receiver" of type GEMMACH (*GEMMACH sets up a pointer to the GEMMACH struct)that then has the method  "read" or whatever.
//So it goes func (receiver *struct) funcname(input inputtype) NextData (which is a function to get the next timestep data from preproc.go) {stuff the function does}.
//This format lets you call the function as w.read(varName) etc.
//Added "noChemHour" - this tells it that there is not an hour and will read-in netCDF files without a time dimension

//SB 20221126: lat and lon stepping through non-uniformly; might be an issue
//not looking into time currently as reading one time at a time: using NOCHEM bool
//level conversions could use a and b parameters defining the levels; check with TR

//TR221118 - We need to flip the gemmach variables as they are read, since they have level 0 = top, and inMAP is level 0 = surface
//Flip is a variable which calls a function that will invert the Z axis of a densearray
//Not entirely sure why it had to be built this way - copied structure from geoschem.go
var flip = func() func(NextData) NextData {
	return func(in NextData) NextData {
		return func() (*sparse.DenseArray, error) {
			data, err := in()
			if err != nil {
				return nil, err
			}
			kmax := data.Shape[0] - 1
			out := sparse.ZerosDense(data.Shape...)
			for k := 0; k < out.Shape[0]; k++ {
				for j := 0; j < out.Shape[1]; j++ {
					for i := 0; i < out.Shape[2]; i++ {
						//Flip it good!
						dataflipped := data.Get(kmax-k, j, i)
						out.Set(dataflipped, k, j, i)
					}
				}
			}
			return out, nil
		}
	}
}

//This function will read the data with the Z axis inverted
//Trying to call the nextDataNCF function outside the return statement did not work - found this work-around in geoschem's preprocessor
func (w *GEMMACH) flipread(varName string) NextData {
	if w.gemnoChemHour {
		flipfunc := flip()
		return flipfunc(nextDataNCF(w.gem_out, gemFormat, varName, w.start, w.end, w.recordDelta, w.fileDelta, readNCFNoHour, w.msgChan))
	} else {
		flipfunc := flip()
		return flipfunc(nextDataNCF(w.gem_out, gemFormat, varName, w.start, w.end, w.recordDelta, w.fileDelta, readNCF, w.msgChan))
	}
}

//This function will read a group of data using the "alt" method with the Z axis inverted
//20230823 - Now reads in invalt() so that the pollutant and alt are on the same level, rather than opposites.
func (w *GEMMACH) flipreadGroupAlt(varGroup map[string]float64) NextData {
	if w.gemnoChemHour {
		flipfunc := flip()
		return flipfunc(nextDataGroupAltNCF(w.gem_out, gemFormat, varGroup, w.invALT(), w.start, w.end, w.recordDelta, w.fileDelta, readNCFNoHour, w.msgChan))
	} else {
		flipfunc := flip()
		return flipfunc(nextDataGroupAltNCF(w.gem_out, gemFormat, varGroup, w.invALT(), w.start, w.end, w.recordDelta, w.fileDelta, readNCF, w.msgChan))
	}
}

//This function reads a group of data with the Z axis inverted
func (w *GEMMACH) flipreadGroup(varGroup map[string]float64) NextData {
	if w.gemnoChemHour {
		flipfunc := flip()
		return flipfunc(nextDataGroupNCF(w.gem_out, gemFormat, varGroup, w.start, w.end, w.recordDelta, w.fileDelta, readNCFNoHour, w.msgChan))
	} else {
		flipfunc := flip()
		return flipfunc(nextDataGroupNCF(w.gem_out, gemFormat, varGroup, w.start, w.end, w.recordDelta, w.fileDelta, readNCF, w.msgChan))
	}
}

func (w *GEMMACH) read(varName string) NextData {
	if w.gemnoChemHour {
		return nextDataNCF(w.gem_out, gemFormat, varName, w.start, w.end, w.recordDelta, w.fileDelta, readNCFNoHour, w.msgChan)
	}
	return nextDataNCF(w.gem_out, gemFormat, varName, w.start, w.end, w.recordDelta, w.fileDelta, readNCF, w.msgChan)
}

func (w *GEMMACH) readGroup(varGroup map[string]float64) NextData {
	if w.gemnoChemHour {
		return nextDataGroupNCF(w.gem_out, gemFormat, varGroup, w.start, w.end, w.recordDelta, w.fileDelta, readNCFNoHour, w.msgChan)
	}
	return nextDataGroupNCF(w.gem_out, gemFormat, varGroup, w.start, w.end, w.recordDelta, w.fileDelta, readNCF, w.msgChan)
}

func (w *GEMMACH) rdpsread(varName string) NextData {
	if w.gemnoChemHour {
		return nextDataNCF(w.gem_rdps, gemFormat, varName, w.start, w.end, w.recordDelta, w.fileDelta, readNCFNoHour, w.msgChan)
	}
	return nextDataNCF(w.gem_rdps, gemFormat, varName, w.start, w.end, w.recordDelta, w.fileDelta, readNCF, w.msgChan)
}

//For these three functions - changed from ALT to AF as GEMMACH
//outputs do not contain ALT, but do contain AF
// Nx helps fulfill the Preprocessor interface by returning
// the number of grid cells in the West-East direction.
//TR20230823 - Trying to reduce calls, going to define global variables Nx,Ny,Nz.
func (w *GEMMACH) Nx() (int, error) {
	f, ff, err := ncfFromTemplate(w.gem_out, gemFormat, w.start)
	if err != nil {
		return -1, fmt.Errorf("nx: %v", err)
	}
	defer f.Close()
	Nx = ff.Header.Lengths("AF")[3]
	return ff.Header.Lengths("AF")[3], nil
}

// Ny helps fulfill the Preprocessor interface by returning
// the number of grid cells in the South-North direction.
func (w *GEMMACH) Ny() (int, error) {
	f, ff, err := ncfFromTemplate(w.gem_out, gemFormat, w.start)
	if err != nil {
		return -1, fmt.Errorf("ny: %v", err)
	}
	defer f.Close()
	Ny = ff.Header.Lengths("AF")[2]
	return ff.Header.Lengths("AF")[2], nil
}

// Nz helps fulfill the Preprocessor interface by returning
// the number of grid cells in the below-above direction.
func (w *GEMMACH) Nz() (int, error) {
	f, ff, err := ncfFromTemplate(w.gem_out, gemFormat, w.start)
	if err != nil {
		return -1, fmt.Errorf("nz: %v", err)
	}
	defer f.Close()
	Nz = ff.Header.Lengths("AF")[1]
	return ff.Header.Lengths("AF")[1], nil
}

// Shp tells us the shape of the output unstaggered array
func (w *GEMMACH) Shp(staggered bool) ([]int, error) {
	f, ff, err := ncfFromTemplate(w.gem_out, gemFormat, w.start)
	if err != nil {
		return []int{-1}, fmt.Errorf("nz: %v", err)
	}
	defer f.Close()
	//Read in dimensions. Hoping this gives performance boost over
	//reading in a variable with the correct shape, but still has I/O calls.
	Nz := ff.Header.Lengths("AF")[1]
	Ny := ff.Header.Lengths("AF")[2]
	Nx := ff.Header.Lengths("AF")[3]
	//Define shape based on Nz,Ny,Nx
	//Add one to Nz if it is staggered
	if staggered {
		shp := []int{Nz + 1, Ny, Nx}
		return shp, err
	} else {
		shp := []int{Nz, Ny, Nx}
		return shp, err
	}

}

// PBLH helps fulfill the Preprocessor interface by returning
// planetary boundary layer height [m].
func (w *GEMMACH) PBLH() NextData { return w.read("H") }

// Height helps fulfill the Preprocessor interface by returning
// layer heights above ground level calculated based on geopotential height.
//TRSB - Used geopotential height directly for height, converted from decametres
// Needs to be on staggered grid (level 4) for further processing
func (w *GEMMACH) Height() NextData {
	hhFunc := w.flipread("GZ") // Geopotential height
	return func() (*sparse.DenseArray, error) {
		HH, err := hhFunc()
		if err != nil {
			return nil, err
		}
		staggered := true
		shp, err := w.Shp(staggered)
		if err != nil {
			return nil, err
		}
		out := sparse.ZerosDense(shp...)
		for k := 0; k < out.Shape[0]; k++ {
			for j := 0; j < out.Shape[1]; j++ {
				for i := 0; i < out.Shape[2]; i++ {
					//Grab the appropriate value (staggered - so 2k)
					HHconv := HH.Get(2*k, j, i) * 10
					if HHconv < 0 {
						HHconv = 1e-10
					}
					out.Set(HHconv, k, j, i)

				}
			}
		}
		//out, err = w.flip(out)
		return out, nil
	}
}

// ALT helps fulfill the Preprocessor interface by returning
// inverse air density [m3/kg].
func (w *GEMMACH) ALT() NextData {
	rhoFunc := w.flipread("RHO") // Density kg/m³
	return func() (*sparse.DenseArray, error) {
		RHO, err := rhoFunc()
		if err != nil {
			return nil, err
		}
		staggered := false
		shp, err := w.Shp(staggered)
		if err != nil {
			return nil, err
		}
		out := sparse.ZerosDense(shp...)
		for k := 0; k < out.Shape[0]; k++ {
			for j := 0; j < out.Shape[1]; j++ {
				for i := 0; i < out.Shape[2]; i++ {
					//On "Level 2" - same as level 1 but with an extra surface layer
					//So, we will take the average of the first two layers then the rest will be offset by 1
					//(since layer 0 is the surface for "level 2" but not "level 1")
					if k == 0 {
						ALTprime := 1. / ((RHO.Get(k, j, i) + RHO.Get(k+1, j, i)) / 2.)
						out.Set(ALTprime, k, j, i)
						continue
					}
					ALTprime := 1. / (RHO.Get(k+1, j, i))
					out.Set(ALTprime, k, j, i)
					continue
				}
			}
		}
		return out, nil
	}
}

//invALT() reads in alt without flipping it, for use in the readgroupalt function.
//Previously, this function was multiplying alt from the bottom level with the pollutant
//concentration at the top.
func (w *GEMMACH) invALT() NextData {
	rhoFunc := w.read("RHO") // Density kg/m³
	return func() (*sparse.DenseArray, error) {
		RHO, err := rhoFunc()
		if err != nil {
			return nil, err
		}
		staggered := false
		shp, err := w.Shp(staggered)
		if err != nil {
			return nil, err
		}
		out := sparse.ZerosDense(shp...)
		for k := 0; k < out.Shape[0]; k++ {
			for j := 0; j < out.Shape[1]; j++ {
				for i := 0; i < out.Shape[2]; i++ {
					//On "Level 2" - same as level 1 but with an extra surface layer
					//So, we will take the average of the first two layers then the rest will be offset by 1
					//(since layer 0 is the surface for "level 2" but not "level 1")
					if k == 0 {
						ALTprime := 1. / ((RHO.Get(k, j, i) + RHO.Get(k+1, j, i)) / 2.)
						out.Set(ALTprime, k, j, i)
						continue
					}
					ALTprime := 1. / (RHO.Get(k+1, j, i))
					out.Set(ALTprime, k, j, i)
					continue
				}
			}
		}
		return out, nil
	}
}

// U helps fulfill the Preprocessor interface by returning U windspeed
//Converting from knots to m/s
// West-East wind speed [m/s].
func (w *GEMMACH) U() NextData {
	uuFunc := w.flipread("UU") // South-North wind speed [kts].
	return func() (*sparse.DenseArray, error) {
		UU, err := uuFunc()
		if err != nil {
			return nil, err
		}
		staggered := false
		shp, err := w.Shp(staggered)
		if err != nil {
			return nil, err
		}
		out := sparse.ZerosDense(shp...)
		for k := 0; k < out.Shape[0]; k++ {
			for j := 0; j < out.Shape[1]; j++ {
				for i := 0; i < out.Shape[2]; i++ {
					//Staggered, so we take the average from cells top and bottom
					//463/900 converts from kts to m/s
					UUprime := ((UU.Get(k, j, i) + UU.Get(k+1, j, i)) / 2) * 463.0 / 900.0
					out.Set(UUprime, k, j, i)
				}
			}
		}
		return out, nil
	}
}

// V helps fulfill the Preprocessor interface by returning
// South-North wind speed [m/s].
func (w *GEMMACH) V() NextData {
	vvFunc := w.flipread("VV") // South-North wind speed [m/s].
	return func() (*sparse.DenseArray, error) {
		VV, err := vvFunc()
		if err != nil {
			return nil, err
		}
		staggered := false
		shp, err := w.Shp(staggered)
		if err != nil {
			return nil, err
		}
		out := sparse.ZerosDense(shp...)
		for k := 0; k < out.Shape[0]; k++ {
			for j := 0; j < out.Shape[1]; j++ {
				for i := 0; i < out.Shape[2]; i++ {
					//Staggered, so we take the average from cells top and bottom
					//463/900 converts from kts to m/s
					VVprime := ((VV.Get(k, j, i) + VV.Get(k+1, j, i)) / 2) * 463.0 / 900.0
					out.Set(VVprime, k, j, i)
				}
			}
		}
		return out, nil
	}
}

func (w *GEMMACH) W() NextData {
	wwFunc := w.flipread("WW") // windspeed  in pa/s
	gzFunc := w.Height()       //Geopotential height (m)
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
		//WW is on level 1 so we need to convert it to a staggered grid for InMap compatability.
		out := sparse.ZerosDense(GZ.Shape...)
		for k := 0; k < out.Shape[0]; k++ {
			for j := 0; j < out.Shape[1]; j++ {
				for i := 0; i < out.Shape[2]; i++ {
					//Calculate the slope of the height/pressure curve
					//For the topmost cell, we will just use the slope from below
					if k >= out.Shape[0]-2 {
						//print(GZ.Get(k, j, i), GZ.Get(k-1, j, i))
						//slope := (GZ.Get(k, j, i) - GZ.Get(k-1, j, i)) /
						//	(PP.Get(k, j, i) - PP.Get(k-1, j, i))
						//Convert as pa/s * m/pa = m/s
						WWprime := out.Get(k-1, j, i)
						out.Set(WWprime, k, j, i)
						continue
					}
					slope := (GZ.Get(k+1, j, i) - GZ.Get(k, j, i)) /
						(PP.Get(k+1, j, i) - PP.Get(k, j, i))
					//Convert as pa/s * m/pa = m/s
					WWprime := WW.Get(k, j, i) * slope
					out.Set(WWprime, k, j, i)
					continue
				}
			}
		}
		return out, nil
	}
}

// AVOC helps fulfill the Preprocessor interface.
func (w *GEMMACH) AVOC() NextData { return w.flipreadGroupAlt(w.aVOC) }

// BVOC helps fulfill the Preprocessor interface.
func (w *GEMMACH) BVOC() NextData { return w.flipreadGroupAlt(w.bVOC) }

// NOx helps fulfill the Preprocessor interface.
func (w *GEMMACH) NOx() NextData { return w.flipreadGroupAlt(w.nox) }

// SOx helps fulfill the Preprocessor interface.
func (w *GEMMACH) SOx() NextData { return w.flipreadGroupAlt(w.sox) }

// NH3 helps fulfill the Preprocessor interface.
func (w *GEMMACH) NH3() NextData { return w.flipreadGroupAlt(w.nh3) }

// ASOA helps fulfill the Preprocessor interface.
// We are going to make it return an error to indicate that it needs
// to be allocated - e.g. that this value is tSOA
func (w *GEMMACH) ASOA() NextData {
	asoaFunc := w.flipreadGroupAlt(w.aSOA)
	return func() (*sparse.DenseArray, error) {
		TSOA, err := asoaFunc()
		if err != nil {
			return nil, err
		}
		//The error message tells the preproc function that this is
		err = fmt.Errorf("tSOA")
		return TSOA, err
	}
}

// BSOA helps fulfill the Preprocessor interface.
//func (w *GEMMACH) BSOA() NextData { return w.flipreadGroupAlt(w.bSOA) }
func (w *GEMMACH) BSOA() NextData {
	bsoaFunc := w.flipreadGroupAlt(w.bSOA)
	return func() (*sparse.DenseArray, error) {
		TSOA, err := bsoaFunc()
		if err != nil {
			return nil, err
		}
		err = fmt.Errorf("tSOA")
		return TSOA, err
	}
}

// PNO helps fulfill the Preprocessor interface.
func (w *GEMMACH) PNO() NextData { return w.flipreadGroupAlt(w.pNO) }

// PS helps fulfill the Preprocessor interface.
func (w *GEMMACH) PS() NextData { return w.flipreadGroupAlt(w.pS) }

// PNH helps fulfill the Preprocessor interface.
//TRSB - changed to "read group" not "read group alt" as value already in ug/kg
//20230811 TR - changed back to ALT as we want it in ug/m3, which is what readgroupalt does.
func (w *GEMMACH) PNH() NextData { return w.flipreadGroupAlt(w.pNH) }

// TotalPM25 helps fulfill the Preprocessor interface.
func (w *GEMMACH) TotalPM25() NextData { return w.flipreadGroup(w.totalPM25) }

// SurfaceHeatFlux helps fulfill the Preprocessor interface
// by returning heat flux at the surface [W/m2].
//FC5 is the aggregate heat flux
func (w *GEMMACH) SurfaceHeatFlux() NextData { return w.read("FC5") }

// UStar helps fulfill the Preprocessor interface
// by returning friction velocity [m/s].
func (w *GEMMACH) UStar() NextData { return w.read("UE") }

// T helps fulfill the Preprocessor interface by
// returning temperature [K].
// TT converted from "level 2" to "level 1", will fail if TT doesn't have 1 extra height layer.
func (w *GEMMACH) T() NextData {
	TTFunc := w.flipread("TT") // Temperature °C
	return func() (*sparse.DenseArray, error) {
		TT, err := TTFunc()
		if err != nil {
			return nil, err
		}
		staggered := false
		shp, err := w.Shp(staggered)
		if err != nil {
			return nil, err
		}
		out := sparse.ZerosDense(shp...)
		for k := 0; k < out.Shape[0]; k++ {
			for j := 0; j < out.Shape[1]; j++ {
				for i := 0; i < out.Shape[2]; i++ {
					//On "Level 2" - same as level 1 but with an extra surface layer
					//So, we will take the average of the first two layers then the rest will be offset by 1
					//(since layer 0 is the surface for "level 2" but not "level 1")
					if k == 0 {
						TTprime := (TT.Get(k, j, i)+TT.Get(k+1, j, i))/2 + 273.15
						out.Set(TTprime, k, j, i)
						continue
					}
					TTprime := (TT.Get(k+1, j, i)) + 273.15
					out.Set(TTprime, k, j, i)
					continue
				}
			}
		}
		return out, nil
	}
}

// P helps fulfill the Preprocessor interface
// by returning pressure [Pa] as pressure level converted
// Here we convert our pressure index to a 3D variable
func (w *GEMMACH) P() NextData {
	ppfunc := w.read("pressure")
	// pnhFunc := w.PNH() // windspeed  in pa/s
	return func() (*sparse.DenseArray, error) {
		PP, err := ppfunc()
		if err != nil {
			return nil, err
		}
		//We will set the shape equal to the U component of windspeed
		// PNH, err := pnhFunc()
		// if err != nil {
		// 	return nil, err
		// }
		// out := sparse.ZerosDense(PNH.Shape...)
		staggered := false
		shp, err := w.Shp(staggered)
		if err != nil {
			return nil, err
		}
		out := sparse.ZerosDense(shp...)
		//slope := (GZ.Get(k+1, j, i) - GZ.Get(k, j, i)) /
		//(PP.Get(k+1, j, i) - PP.Get(k, j, i))
		//This is probably slower than it needs to be - going element by element.
		//Setting the value at height K for each element. Flipping here so that surface at 0
		kmax := PP.Shape[0] - 1
		for k := 0; k < out.Shape[0]; k++ {
			for j := 0; j < out.Shape[1]; j++ {
				for i := 0; i < out.Shape[2]; i++ {
					out.Set(PP.Get(kmax-k), k, j, i)
				}
			}
		}

		return out.ScaleCopy(100.0), nil
	}
}

// HO helps fulfill the Preprocessor interface
// by returning hydroxyl radical concentration [ppmv].
func (w *GEMMACH) HO() NextData { return w.flipreadGroup(w.ho) }

// H2O2 helps fulfill the Preprocessor interface
// by returning hydrogen peroxide concentration [ppmv].
func (w *GEMMACH) H2O2() NextData { return w.flipreadGroup(w.h2o2) }

//SB 20221127
//Start Copy-Adapt-VF-v0
//creating functions to take in values needed for running landuse dimensions from Gem_geophy
// xCenters returns the x-coordinates of the grid points.

//SB 20221127: Copy-Adapt-VF-v0
//gem_geophy holds information about VF on the domain.
//could bring in each VF using this
type gem_geophy struct {
	data   *rtree.Rtree
	dx, dy float64
	xo, yo float64
	nx, ny int
}

//SB 20221127: category int being used to store known VF in each cell
// adding vf as a measure of vf fraction in each cell
type gemGridCell struct {
	geom.Polygon
	vf       float64
	category int
}

// DX returns the longitude grid spacing.
func DX(file *cdf.File) (float64, error) {
	//return float64(file.Header.GetAttribute("", "delta_rlon").([]float32)[0]), nil
	return float64(file.Header.GetAttribute("", "delta_x").([]float64)[0]), nil
	//return w.chemAttribute("delta_rlon")
}

// DY returns the latitude grid spacing.
func DY(file *cdf.File) (float64, error) {
	//return float64(file.Header.GetAttribute("", "delta_rlat").([]float32)[0]), nil
	return float64(file.Header.GetAttribute("", "delta_y").([]float64)[0]), nil
	//return w.chemAttribute("delta_rlat")
}

func xCenters(file *cdf.File) ([]float64, error) {
	//something using the interface Reader
	//var v = "rLON_var" "x_var"
	var v = "x"
	dims := file.Header.Lengths(v)
	// //ny := dims[0]
	//nx := dims[1]
	nx := dims[0]
	data, err := readNCFNoHour(v, file, 0)
	if err != nil {
		// If variable not in file, try all lowercase.
		data, err = readNCFNoHour(strings.ToLower(v), file, 0)
		if err != nil {
			return nil, err
		}
	}
	return data.Elements[:nx], nil
}

// yCenters returns the y-coordinates of the grid points.
func yCenters(file *cdf.File) ([]float64, error) {
	//something using the interface Reader
	//var v = "rLAT_var"
	var v = "y"
	dims := file.Header.Lengths(v)
	ny := dims[0]
	//nx := dims[1]
	//yCenter := make([]float64, ny)
	data, err := readNCFNoHour(v, file, 0)
	if err != nil {
		// If variable not in file, try all lowercase.
		data, err = readNCFNoHour(strings.ToLower(v), file, 0)
		if err != nil {
			return nil, err
		}
	}
	return data.Elements[:ny], nil
}

// SB 20221127:
// readGem_geophy reads data from an Gem_geophy file:
// It may be downlaodable from:
// https://hpfx.collab.science.gc.ca/~smco810/Deliveries/collab/rcm/D001/
// Gem_geophy also needed input variables rLAT_var rLON_var and delta_rlat delta_rlon
// created those variables in r and python and inserted in netcdf file
// The file must be converted from netcdf version 4 to version 3 before
// use by this function.
// This can be done using the command:
// nccopy -k classic Gem_geophy.nc Gem_geophy1.nc
// name change is necessary to ensure code works
// added a new input var i of type int to mark vegetation type from nc file
func (w *GEMMACH) readgem_geophy(file *cdf.File, i string) (*gem_geophy, error) {

	//taking in global attributes from gem_geophy
	dx := float64(file.Header.GetAttribute("", "delta_x").([]float64)[0])
	//dxStr := file.Header.GetAttribute("", "delta_rlon").(string)
	dy := float64(file.Header.GetAttribute("", "delta_y").([]float64)[0])

	//originally olson land map was running here
	//it was taking in an integer representation of largest land use in very refined grid cell
	//INMAP uses coarser grids that do not exactly align with the underlying input grid
	//so needs to calculate the largest land use fraction in the new grid
	//code below only takes in the integer representations and stores them in the original grid
	dims := file.Header.Lengths("VF" + fmt.Sprintf("%v", i))
	ny := dims[0]
	nx := dims[1]

	r := file.Reader("VF"+fmt.Sprintf("%v", i), nil, nil)

	buf := r.Zero(-1).([]float32)
	o := &gem_geophy{
		data: rtree.NewTree(25, 50),
		xo:   w.xCenters[0],
		yo:   w.yCenters[0],
		dx:   dx,
		dy:   dy,
		nx:   nx,
		ny:   ny,
	}

	if _, err := r.Read(buf); err != nil {
		return nil, fmt.Errorf("inmap: reading Gem_geophy: %v", err)
	}

	for iy := 0; iy < o.ny; iy++ {
		for ix := 0; ix < o.nx; ix++ {
			//dy := w.dy
			if iy == 0 {
				dy = ((w.yCenters[iy+1] - w.yCenters[iy]) - w.dy/2) * 2
			}
			if iy == len(w.yCenters)-1 {
				dy = ((w.yCenters[iy] - w.yCenters[iy-1]) - w.dy/2) * 2
			}
			//for i, x := range w.xCenters {
			//dx := w.dx
			if ix == 0 {
				dx = ((w.xCenters[ix+1] - w.xCenters[ix]) - w.dx/2) * 2
			}
			if ix == len(w.xCenters)-1 {
				dx = ((w.xCenters[ix] - w.xCenters[ix-1]) - w.dx/2) * 2
			}
			x0 := w.xCenters[ix] - dx/2
			x1 := w.xCenters[ix] + dx/2
			y0 := w.yCenters[iy] - dy/2
			y1 := w.yCenters[iy] + dy/2
			//x0 := o.xo + o.dx*float64(ix)
			//x1 := o.xo + o.dx*float64(ix+1)
			//y0 := o.yo + o.dy*float64(iy)
			//y1 := o.yo + o.dy*float64(iy+1)
			v, _ := strconv.Atoi(i)
			c := gemGridCell{
				Polygon: geom.Polygon{{
					{X: x0, Y: y0},
					{X: x1, Y: y0},
					{X: x1, Y: y1},
					{X: x0, Y: y1},
				}},
				category: int(v),
				vf:       float64(buf[iy*nx+ix]),
			}
			o.data.Insert(c)
		}
	}
	return o, nil
}

// fractions returns the fraction of land use types within the given polygon.
//o is the landuse data - gets put into 'c' to go through the 26 VF types. p is the geometry we are putting it on.
func fractions(p geom.Polygon, o []gem_geophy) map[int]float64 {
	out := make(map[int]float64)
	for i := 0; i < 26; i++ {
		for _, cI := range o[i].data.SearchIntersect(p.Bounds()) {
			c := cI.(gemGridCell)
			isect := p.Intersection(c)
			//print(c.Area())
			if isect != nil { // Add the area of the intersection between the VF grid and the INMAP grid for the given VF type
				//out[c.category] += isect.Area()
				out[c.category] += isect.Area() * c.vf //TR20221130 - function wasn't accounting for actual veg frac, needed to multiply.
			}
		}
	}
	//TR20221130 - moved this out of the above for loop - we want to divide the area of each fraction by the total area.
	a := p.Area()

	for cat := range out {
		out[cat] /= a //Divide the area for each category by the total polygon area
	}

	return out
}

// largestLandUse returns the land use index with the largest area
// in each grid cell when given a gem_geophy file.
func (w *GEMMACH) largestLandUse(file *cdf.File) (*sparse.DenseArray, error) {
	var VFs [26]string
	var o = make([]gem_geophy, 26)

	//var p gem_geophy
	//var err error
	for i := 0; i < 26; i++ {
		if i < 9 {
			VFs[i] = "0" + fmt.Sprintf("%v", i+1)
		} else {
			VFs[i] = fmt.Sprintf("%v", i+1)
		}
		p, err := w.readgem_geophy(file, VFs[i])
		if err != nil {
			return nil, err
		}
		//o = append(o, *p)
		o[i] = *p
		//fmt.Println(o[i])
		//i++
	}

	//new grid being created and fractions of VF being calculated here
	out := sparse.ZerosDense(len(w.yCenters), len(w.xCenters))
	for j, y := range w.yCenters {
		dy := w.dy
		if j == 0 {
			dy = ((w.yCenters[j+1] - y) - w.dy/2) * 2
		}
		if j == len(w.yCenters)-1 {
			dy = ((y - w.yCenters[j-1]) - w.dy/2) * 2
		}
		for i, x := range w.xCenters {
			dx := w.dx
			if i == 0 {
				dx = ((w.xCenters[i+1] - x) - w.dx/2) * 2
			}
			if i == len(w.xCenters)-1 {
				dx = ((x - w.xCenters[i-1]) - w.dx/2) * 2
			}
			x0 := x - dx/2
			x1 := x + dx/2
			y0 := y - dy/2
			y1 := y + dy/2
			p := geom.Polygon{{
				{X: x0, Y: y0},
				{X: x1, Y: y0},
				{X: x1, Y: y1},
				{X: x0, Y: y1},
				//{X: x0, Y: y0},
			}}

			fractions := fractions(p, o)
			if len(fractions) == 0 {
				return nil, fmt.Errorf("no land use information available for polygon %+v", p)
			}

			maxCat := math.MinInt32
			maxVal := math.Inf(-1)
			for c, v := range fractions {
				if v > maxVal {
					maxVal = v
					maxCat = c
				}
			}

			out.Set(float64(maxCat), j, i)
		}
	}
	return out, nil
}

// SeinfeldLandUse helps fulfill the Preprocessor interface
// by returning land use categories as
// specified in github.com/ctessum/atmos/seinfeld.
func (w *GEMMACH) SeinfeldLandUse() NextData {
	//SB 20221128: there is no default index to read unlike WRF
	//like GEOS-Chem need to read in the sparse densearray created for VF types
	//luFunc := w.read("LU_INDEX") // USGS land use index
	return w.GEMSeinfeldLandUse(w.landUse)
}

func (w *GEMMACH) GEMSeinfeldLandUse(landUse *sparse.DenseArray) NextData {
	pblhFunc := w.PBLH()
	//SB: no need to use flip because no Z component
	//fnFunc := w.flipread(landuse)
	return func() (*sparse.DenseArray, error) {
		//SB 20221128: unlike GEOS-Chem, no separate snowfrac
		//snowFrac, err := snowFunc() // Fraction land covered by snow
		// if err != nil {
		// 	return nil, err
		// }
		//We will set the shape equal to the U component of windspeed
		PBLH, err := pblhFunc()
		if err != nil {
			return nil, err
		}
		o := sparse.ZerosDense(PBLH.Shape...)
		//o := sparse.ZerosDense(snowFrac.Shape...)
		for j := 0; j < o.Shape[0]; j++ {
			for i := 0; i < o.Shape[1]; i++ {
				//Convert GEM landuse to Seinfeld. Need to subtract 1 as index[0] = VF1
				oval := GEMseinfeld[f2i(landUse.Get(j, i))-1]
				//print(float64(oval))
				o.Set(float64(oval), j, i)
			}
		}
		return o, nil
	}
}

// func GEMSeinfeldLandUse(luFunc NextData) NextData {
// 	return func() (*sparse.DenseArray, error) {
// 		lu, err := luFunc() // USGS land use index
// 		if err != nil {
// 			return nil, err
// 		}
// 		o := sparse.ZerosDense(lu.Shape...)
// 		for j := 0; j < lu.Shape[0]; j++ {
// 			for i := 0; i < lu.Shape[1]; i++ {
// 				o.Set(float64(GEMseinfeld[f2i(lu.Get(j, i))]), j, i)
// 			}
// 		}
// 		return o, nil
// 	}
// }

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
	//SB 20221128: there is no default index to read unlike WRF
	//like GEOS-Chem need to read in the sparse densearray created for VF types
	//luFunc := w.read("LU_INDEX") // USGS land use index
	return w.GEMweselyLandUse(w.landUse)
}

// func (w *GEMMACH) WeselyLandUse() NextData {
// 	luFunc := w.read("VF_INDEX") // USGS land use index
// 	return GEMWeselyLandUse(luFunc)
// }

func (w *GEMMACH) GEMweselyLandUse(landUse *sparse.DenseArray) NextData {
	pblhFunc := w.PBLH() //TR20221130 changed to PBLH as it is a 2D variable
	return func() (*sparse.DenseArray, error) {
		// snowFrac, err := snowFunc() // Fraction land covered by snow
		// if err != nil {
		// 	return nil, err
		// }
		PBLH, err := pblhFunc()
		if err != nil {
			return nil, err
		}
		o := sparse.ZerosDense(PBLH.Shape...)
		for j := 0; j < o.Shape[0]; j++ {
			for i := 0; i < o.Shape[1]; i++ {
				//TR 20221130 - Needed to subtract 1 from the landuse value and changed to 2D variable.
				oval := GEMwesely[f2i(landUse.Get(j, i))-1]
				o.Set(float64(oval), j, i)
			}
		}
		return o, nil
	}
}

// func (w *GEMMACH) GEMWeselyLandUse(luFunc NextData) NextData {
// 	return func() (*sparse.DenseArray, error) {
// 		lu, err := luFunc() // USGS land use index
// 		if err != nil {
// 			return nil, err
// 		}
// 		o := sparse.ZerosDense(lu.Shape...)
// 		for j := 0; j < lu.Shape[0]; j++ {
// 			for i := 0; i < lu.Shape[1]; i++ {
// 				o.Set(float64(GEMVFwesely[f2i(lu.Get(j, i))]), j, i)
// 			}
// 		}
// 		return o, nil
// 	}
// }

// GEMVFwesely lookup table to go from GEM vegetation fraction to land classes for
// gas dry deposition.
//TR - need to convert from GEMMACH wesely to INMAP wesely
var GEMwesely = []wesely1989.LandUseCategory{
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
	//LUIndexFunc := w.read("VF_INDEX") //GEM VF land use index
	return w.GEMZ0(w.landUse)
}

func (w *GEMMACH) GEMZ0(landUse *sparse.DenseArray) NextData {
	pblhFunc := w.PBLH()
	return func() (*sparse.DenseArray, error) {
		PBLH, err := pblhFunc()
		if err != nil {
			return nil, err
		}
		// luIndex, err := LUIndexFunc()
		// if err != nil {
		// 	return nil, err
		// }
		zo := sparse.ZerosDense(PBLH.Shape...)
		for i, lu := range landUse.Elements {
			zo.Elements[i] = GEMz0[f2i(lu)-1] // roughness length [m]
			//oval := GEMwesely[f2i(landUse.Get(j, i))-1]
			//	o.Set(float64(oval), j, i)

		}
		return zo, nil
	}
}

// GEMz0 holds Roughness lengths for GEM land classes ([m]), from WRF file
// TR - these are the GEMMACH definitions
var GEMz0 = []float64{0.001, 0.0003, 0.001, 1.5, 3.5, 1,
	2.0, 3.0, .80, .05, .15, .15, .02, .08, 0.08, .08, .35,
	.25, 0.1, .08, 1.35, .01, .05, .05, 1.5, .05}

// QRain helps fulfill the Preprocessor interface by
// returning rain mass fraction.
//TR -Sahil to do
//func (w *GEMMACH) QRain() NextData { return w.read("QRAIN") }

// QRain helps fulfill the Preprocessor interface by returning
// rain mass fraction based on the GEM precipitation rate [m s-1]
// and the assumption (from the EMEP model wet deposition algorithm)
// that raindrops are falling at 5 m/s.
func (w *GEMMACH) QRain() NextData {
	PrecipFunc := w.read("RT") // Total precipitation rate (m/s)
	altFunc := w.ALT()         // 1/density
	cloudFunc := w.CloudFrac() //Fraction of cloud in each grid unit.
	return func() (*sparse.DenseArray, error) {
		precipgem, err := PrecipFunc()
		if err != nil {
			return nil, err
		}
		alt, err := altFunc()
		if err != nil {
			return nil, err
		}
		cloud, err := cloudFunc()
		if err != nil {
			return nil, err
		}
		const Vdr = 5.0 // droplet velocity [m/s] EMEP wet dep algorithm assumption
		qRain := sparse.ZerosDense(alt.Shape...)
		// Need to define a 3D precipitation variable.
		for j := 0; j < qRain.Shape[1]; j++ {
			for i := 0; i < qRain.Shape[2]; i++ {
				// Allocate rainfall in every grid cell with rainfall from the surface to
				// the highest cloud, or the top of the column if no clouds present.
				maxk := float64(qRain.Shape[0]) + 1.
				for k := 0; k < qRain.Shape[0]; k++ {
					if cloud.Get(k, j, i) != 0 {
						maxk = float64(k)
					}
				}
				for k := 0; k < qRain.Shape[0]; k++ {
					//Allocate evenly to all cells lower than maxk (since maxk is +1, don't do it)
					precip3D := 0.
					if float64(k) < maxk {
						precip3D = precipgem.Get(j, i) / maxk
					} // From EMEP algorithm: P = QRAIN * Vdr * ρgas => QRAIN = P / Vdr / ρgas
					// [kg m-2 s-1] / [m s-1] * [m3 kg-1]
					q := precip3D / Vdr * alt.Get(k, j, i)
					qRain.Set(q, k, j, i)
				}

			}
		}

		return qRain, nil
	}
}

// CloudFrac helps fulfill the Preprocessor interface
// by returning the fraction of each grid cell filled
// with clouds [volume/volume].
func (w *GEMMACH) CloudFrac() NextData {
	fnFunc := w.flipread("FN") // Cloud Fraction (V/V)
	// pnhFunc := w.PNH()         //Importing for shape. Probably a better way to do this.
	return func() (*sparse.DenseArray, error) {
		FN, err := fnFunc()
		if err != nil {
			return nil, err
		}
		// PNH, err := pnhFunc()
		// if err != nil {
		// 	return nil, err
		// }
		// out := sparse.ZerosDense(PNH.Shape...)
		staggered := false
		shp, err := w.Shp(staggered)
		if err != nil {
			return nil, err
		}
		out := sparse.ZerosDense(shp...)
		for k := 0; k < out.Shape[0]; k++ {
			for j := 0; j < out.Shape[1]; j++ {
				for i := 0; i < out.Shape[2]; i++ {
					//On "Level 2" - same as level 1 but with an extra surface layer
					//So, we will take the average of the first two layers then the rest will be offset by 1
					//(since layer 0 is the surface for "level 2" but not "level 1")
					if k == 0 {
						FNprime := (FN.Get(k, j, i) + FN.Get(k+1, j, i)) / 2
						out.Set(FNprime, k, j, i)
						continue
					}
					FNprime := (FN.Get(k+1, j, i))
					out.Set(FNprime, k, j, i)
					continue
				}
			}
		}
		return out, nil
	}
}

// QCloud helps fulfill the Preprocessor interface by returning
// the mass fraction of cloud water in each grid cell [mass/mass].
func (w *GEMMACH) QCloud() NextData { return w.rdpsread("QC") }

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
