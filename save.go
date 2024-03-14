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
	"encoding/gob"
	"fmt"
	"io"
	"sort"

	"github.com/ctessum/geom"

	//TR20230705 - added msgpack and gzip compression for larger file size compatability
	"compress/gzip"
	"github.com/vmihailenco/msgpack/v5"
	//"math"
)

func init() {
	gob.Register(geom.Polygon{})
	gob.Register(geom.MultiPolygon{})
	gob.Register(&geom.Bounds{})
	gob.Register(geom.Point{})
	gob.Register(geom.MultiPoint{})
	gob.Register(geom.LineString{})
	gob.Register(geom.MultiLineString{})
}

type versionCells struct {
	// DataVersion holds the variable grid data version of the software
	// that saved this data, if any, and should match the VarGridDataVersion
	// global variable.
	DataVersion string
	Cells       []*Cell
}

// Save returns a function that saves the data in d to a gob file
// (format description at https://golang.org/pkg/encoding/gob/).
func Save(w io.Writer) DomainManipulator {
	return func(d *InMAP) error {

		if d.cells.len() == 0 {
			return fmt.Errorf("inmap.InMAP.Save: no grid cells to save")
		}

		// Set the data version so it can be checked when the data is loaded.
		data := versionCells{
			DataVersion: VarGridDataVersion,
			Cells:       d.cells.array(),
		}

		e := gob.NewEncoder(w)

		if err := e.Encode(data); err != nil {
			return fmt.Errorf("inmap.InMAP.Save: %v", err)
		}
		return nil
	}
}

func SaveVarGrid(w io.Writer, d *InMAP) DomainManipulator {
	return func(d *InMAP) error {

		if d.cells.len() == 0 {
			return fmt.Errorf("inmap.InMAP.Save: no grid cells to save")
		}

		// Set the data version so it can be checked when the data is loaded.
		data := versionCells{
			DataVersion: VarGridDataVersion,
			Cells:       d.cells.array(),
		}

		e := gob.NewEncoder(w)

		if err := e.Encode(data); err != nil {
			return fmt.Errorf("inmap.InMAP.Save: %v", err)
		}
		return nil
	}
}

// Load returns a function that loads the data from a previously Saved file
// into an InMAP object.
// func Load(r io.Reader, config *VarGridConfig, emis *Emissions, m Mechanism) DomainManipulator {
// 	return func(d *InMAP) error {
// 		dec := gob.NewDecoder(r)
// 		var data versionCells
// 		if err := dec.Decode(&data); err != nil {
// 			return fmt.Errorf("inmap.InMAP.Load: %v", err)
// 		}
// 		if err := d.initFromCells(data.Cells, emis, config, m); err != nil {
// 			return err
// 		}
// 		if data.DataVersion != VarGridDataVersion {
// 			return fmt.Errorf("InMAP variable grid data version %s is not compatible with "+
// 				"the required version %s", data.DataVersion, VarGridDataVersion)
// 		}
// 		return nil
// 	}
// }
func Load(r io.Reader, config *VarGridConfig, emis *Emissions, m Mechanism) DomainManipulator {
	return func(d *InMAP) error {
		// Try decoding as msgpack
		decMsgpack := msgpack.NewDecoder(r)
		var data versionCells
		if err := decMsgpack.Decode(&data); err == nil {
			if err := d.initFromCells(data.Cells, emis, config, m); err != nil {
				return err
			}
			if data.DataVersion != VarGridDataVersion {
				return fmt.Errorf("InMAP variable grid data version %s is not compatible with the required version %s", data.DataVersion, VarGridDataVersion)
			}
			return nil
		}

		// If msgpack decoding fails, try decoding as gob
		decGob := gob.NewDecoder(r)
		var gobData versionCells
		if err := decGob.Decode(&gobData); err != nil {
			return fmt.Errorf("inmap.InMAP.Load: %v", err)
		}

		if err := d.initFromCells(gobData.Cells, emis, config, m); err != nil {
			return err
		}
		if gobData.DataVersion != VarGridDataVersion {
			return fmt.Errorf("InMAP variable grid data version %s is not compatible with the required version %s", gobData.DataVersion, VarGridDataVersion)
		}
		return nil
	}
}

//SaveBig returns a function that saves the data in d to a JSON file
//with optional gzip compression. Made with chat GPT assistance
func SaveBig(w io.Writer, compress bool) DomainManipulator {
	return func(d *InMAP) error {
		if d.cells.len() == 0 {
			return fmt.Errorf("inmap.InMAP.Save: no grid cells to save")
		}

		var writer io.Writer = w
		if compress {
			gzipWriter := gzip.NewWriter(w)
			defer gzipWriter.Close()
			writer = gzipWriter
		}

		// Set the data version so it can be checked when the data is loaded.
		data := versionCells{
			DataVersion: VarGridDataVersion,
			Cells:       d.cells.array(),
		}

		enc := msgpack.NewEncoder(writer)
		//enc.EncodeOptions = enc.EncodeOptions.WithMarshaler((*CustomFloat64)(nil), (*json.Marshaler)(nil))
		if err := enc.Encode(data); err != nil {
			return fmt.Errorf("inmap.InMAP.SaveBig: %v", err)
		}
		return nil
	}
}

//LoadBig returns a function that loads the data from a previously saved JSON file
//with optional gzip compression into an InMAP object.  Made with chat GPT assistance
// func LoadBig(r io.Reader, config *VarGridConfig, emis *Emissions, m Mechanism) DomainManipulator {
// 	return func(d *InMAP) error {
// 		var reader io.Reader = r
// 		gzipReader, err := gzip.NewReader(r)
// 		if err == nil {
// 			defer gzipReader.Close()
// 			reader = gzipReader
// 		}

// 		// Try decoding with msgpack
// 		dec := msgpack.NewDecoder(reader)
// 		//dec.DecodeOptions = dec.DecodeOptions.WithUnmarshaler((*CustomFloat64)(nil), (*json.Unmarshaler)(nil))
// 		var cells []*Cell
// 		if err := dec.Decode(&cells); err == nil {
// 			// msgpack decoding successful
// 			if err := d.initFromCells(cells, emis, config, m); err != nil {
// 				return err
// 			}
// 			return nil
// 		}

// 		// Reset the reader to the original state
// 		if _, err := r.(io.Seeker).Seek(0, io.SeekStart); err != nil {
// 			return fmt.Errorf("inmap.InMAP.LoadBig: failed to reset reader: %v", err)
// 		}

// 		// Call Load function directly if this fails
// 		return Load(r, config, emis, m)(d)
// 	}
// }

func (d *InMAP) initFromCells(cells []*Cell, emis *Emissions, config *VarGridConfig, m Mechanism) error {
	d.init()
	// Create a list of array indices for each population type.
	d.PopIndices = make(map[string]int)
	for i, p := range config.CensusPopColumns {
		d.PopIndices[p] = i
	}
	d.mortIndices = make(map[string]int)
	mortRateColumns := make([]string, len(config.MortalityRateColumns))
	i := 0
	for m := range config.MortalityRateColumns {
		mortRateColumns[i] = m
		i++
	}
	sort.Strings(mortRateColumns)
	for i, m := range mortRateColumns {
		d.mortIndices[m] = i
	}
	for _, c := range cells {
		d.InsertCell(c, m)
	}

	// Add emissions to new cells.
	// This needs to be called after setNeighbors.
	if err := d.SetEmissionsFlux(emis, m); err != nil {
		return err
	}
	return nil
}

// SaveBig returns a function that saves the data in d to a JSON file with optional gzip compression.
// func SaveBig(w io.Writer, compress bool) DomainManipulator {
// 	return func(d *InMAP) error {
// 		if d.cells.len() == 0 {
// 			return fmt.Errorf("inmap.InMAP.SaveBig: no grid cells to save")
// 		}

// 		var writer io.Writer = w
// 		if compress {
// 			gzipWriter := gzip.NewWriter(w)
// 			defer gzipWriter.Close()
// 			writer = gzipWriter
// 		}

// 		enc := json.NewEncoder(writer)
// 		enc.SetEscapeHTML(false)

// 		// Convert NaN values to a special string representation during encoding
// 		for _, cell := range d.cells.array() {
// 			value := cell.Value.(float64)
// 			if math.IsNaN(value) {
// 				cell.Value = "NaN"
// 			}
// 		}

// 		// Encode using json.NewEncoder
// 		if err := enc.Encode(d.cells.array()); err != nil {
// 			return fmt.Errorf("inmap.InMAP.SaveBig: %v", err)
// 		}
// 		return nil
// 	}
// }

// // LoadBig returns a function that loads the data from a previously saved JSON file
// // with optional gzip compression into an InMAP object.
// func LoadBig(r io.Reader, config *VarGridConfig, emis *Emissions, m Mechanism) DomainManipulator {
// 	return func(d *InMAP) error {
// 		var reader io.Reader = r
// 		gzipReader, err := gzip.NewReader(r)
// 		if err == nil {
// 			defer gzipReader.Close()
// 			reader = gzipReader
// 		}

// 		dec := json.NewDecoder(reader)
// 		dec.UseNumber()

// 		var cells []*Cell
// 		if err := dec.Decode(&cells); err != nil {
// 			// JSON decoding failed, try calling Load function directly
// 			if loadErr := Load(r, config, emis, m)(d); loadErr != nil {
// 				return fmt.Errorf("inmap.InMAP.LoadBig: failed to decode JSON and GOB: %v", err)
// 			}
// 			return nil
// 		}

// 		// Convert special string representation back to NaN values during decoding
// 		for _, cell := range cells {
// 			if strValue, ok := cell.Value.(string); ok && strValue == "NaN" {
// 				cell.Value = math.NaN()
// 			}
// 		}

// 		if err := d.initFromCells(cells, emis, config, m); err != nil {
// 			return err
// 		}
// 		return nil
// 	}
// }

// // // CustomFloat64 is a custom float64 type that handles NaN values during JSON encoding/decoding.
// // type CustomFloat64 float64

// // // MarshalJSON marshals the CustomFloat64 value to JSON, handling NaN values.
// // // Made with chat GPT assistance
// // func (f CustomFloat64) MarshalJSON() ([]byte, error) {
// // 	if math.IsNaN(float64(f)) {
// // 		return []byte("NaN"), nil
// // 	}
// // 	return json.Marshal(float64(f))
// // }

// // // UnmarshalJSON unmarshals the JSON value into the CustomFloat64 type, handling NaN values.
// // // Made with chat GPT assistance
// // func (f *CustomFloat64) UnmarshalJSON(data []byte) error {
// // 	if string(data) == "NaN" {
// // 		*f = CustomFloat64(math.NaN())
// // 		return nil
// // 	}
// // 	var value float64
// // 	err := json.Unmarshal(data, &value)
// // 	if err != nil {
// // 		return err
// // 	}
// // 	*f = CustomFloat64(value)
// // 	return nil
// // }