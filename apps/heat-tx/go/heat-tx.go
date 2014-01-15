// Copyright (c) 2014, Los Alamos National Security, LLC All rights reserved.
//
// This software was produced under U.S. Government contract DE-AC52-06NA25396
// for Los Alamos National Laboratory (LANL), which is operated by Los Alamos
// National Security, LLC for the U.S. Department of Energy. The U.S. Government
// has rights to use, reproduce, and distribute this software.  NEITHER THE
// GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS
// OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If
// software is modified to produce derivative works, such modified software
// should be clearly marked, so as not to confuse it with the version available
// from LANL.
//
// Additionally, redistribution and use in source and binary forms, with or
// without modification, are permitted provided that the following conditions
// are met:
//
// . Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// . Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// . Neither the name of Los Alamos National Security, LLC, Los Alamos National
//   Laboratory, LANL, the U.S. Government, nor the names of its contributors
//   may be used to endorse or promote products derived from this software
//   without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
// CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT
// NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
// SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// LA-CC 10-123

// A simple 2D heat transfer simulation in Go by Samuel K. Gutierrez

// To Build:
// go build heat-tx.go

// To Profile:
// Build
// ./heat-tx -cpuprofile=heat-tx.prof
// go tool pprof ./heat-tx ./heat-tx.prof

// To Plot (gnuplot):
// plot './heat-img.dat' matrix with image

package main

import (
    "flag"
    "strconv"
    "fmt"
    "bufio"
    "os"
    "log"
    "runtime/pprof"
)

// Application constants
const (
    // Application Name
    AppName string = "go-heat-tx"
    // Application version string
    AppVerStr string = "0.1"
    // Max simulation time
    TMax uint64 = 1024
    // Mesh size (x and y)
    N uint64 = 512
    // Thermal conductivity
    ThermCond float64 = 0.6
    // Some constant
    K float64 = 0.4
)

// 2D mesh
type Mesh struct {
    nx, ny uint64 // mesh size in x and y
    cells [][]float64 // mesh cells
}

type SimParams struct {
    // Thermal conductivity
    c float64
    deltaS float64
    // time interval
    deltaT float64
    // Max sim iterations
    tMax uint64
}

type HeatTxSim struct {
    // The meshes
    newMesh, oldMesh *Mesh
    // Simulation parameters
    params *SimParams
}

// NewMesh returns an empty mesh of the specified width and height
func NewMesh(x, y uint64) *Mesh {
    cells := make([][]float64, x)
    for i := range cells {
        cells[i] = make([]float64, y)
    }
    // **remember** unlike C, we can return the address of a local variable.
    // in fact, this returns a fresh instance each time the following code is
    // evaluated - w00t.
    return &Mesh{nx: x, ny: y, cells: cells}
}

// NewSimParams returns a new set of initialized simulation parameters based on
// the provided input.
func NewSimParams(nx uint64, thermCond float64, tMax uint64) *SimParams {
    fmt.Println("o initializing simulation parameters...")
    ds := 1.0 / (float64(nx) + 1.0)
    dt := (ds * ds) / (4.0 * thermCond)
    sp := &SimParams{c: thermCond, deltaS: ds, deltaT: dt, tMax: tMax}
    fmt.Print(sp)
    return sp
}

// Nice SimParams printing
func (p *SimParams) String() string {
    pStr := ""
    pStr += fmt.Sprintf(". max_t: %d\n", p.tMax)
    pStr += fmt.Sprintf(". c: %f\n", p.c)
    pStr += fmt.Sprintf(". delta_s: %f\n", p.deltaS)
    pStr += fmt.Sprintf(". delta_t: %f\n", p.deltaT)
    pStr += "\n"
    return pStr
}

// Nice Mesh String representation
func (m *Mesh) String() string {
    mStr := ""
    for i := range m.cells {
        for j := range m.cells[i] {
            mStr += strconv.FormatFloat(m.cells[i][j], 'e', 1, 64) + " "
        }
        mStr += "\n"
    }
    return mStr
}

func NewHeatTxSim(x, y uint64,  thermCond float64, tMax uint64) *HeatTxSim {
    return &HeatTxSim{params: NewSimParams(x, thermCond, tMax),
                      newMesh: NewMesh(x, y), oldMesh: NewMesh(x, y)}
}

func (m *Mesh) SetInitConds() {
    x0 := m.nx / 2
    y0 := m.ny / 2
    x  := m.nx / 4
    y  := uint64(0)
    radiusErr := int64(1 - x)

    for x >= y {
        m.cells[ x + x0][ y + y0] = K
        m.cells[ x + x0][ y + y0] = K * .50
        m.cells[ y + x0][ x + y0] = K * .60
        m.cells[-x + x0][ y + y0] = K * .70
        m.cells[-y + x0][ x + y0] = K * .80
        m.cells[-x + x0][-y + y0] = K * .70
        m.cells[-y + x0][-x + y0] = K * .60
        m.cells[ x + x0][-y + y0] = K * .50
        m.cells[ y + x0][-x + y0] = K
        y++
        if radiusErr < 0 {
            radiusErr += int64(2 * y + 1)
        } else {
            x--
            radiusErr += int64(2 * (y - x + 1))
        }
    }
}

// Runs the simulation
func (s *HeatTxSim) Run() {
    nx := len(s.oldMesh.cells) - 1
    ny := len(s.oldMesh.cells[0]) - 1
    ds2 := s.params.deltaS * s.params.deltaS
    cdtods2 := (s.params.c * s.params.deltaT) / ds2
    tMax := s.params.tMax;
    newMesh := s.newMesh
    oldMesh := s.oldMesh

    fmt.Println("o starting simulation...")
    for t := uint64(0); t < tMax; t++ {
        if t % 100 == 0 {
            fmt.Println(". starting iteration", t, "of", tMax)
        }
        for i := 1; i < nx; i++ {
            for j := 1; j < ny; j++ {
                newMesh.cells[i][j] =
                    oldMesh.cells[i][j] +
                    (cdtods2 *
                     (oldMesh.cells[i + 1][j] +
                      oldMesh.cells[i - 1][j] -
                      4.0 * oldMesh.cells[i][j] +
                      oldMesh.cells[i][j + 1] +
                      oldMesh.cells[i][j - 1]))
            }
        }
        // swap old and new - this is just a pointer swap
        oldMesh, newMesh = newMesh, oldMesh
        // Constant heat source
        oldMesh.SetInitConds()
    }
}

// Dumps to a text file with the current newMesh cell values
func (s *HeatTxSim) Dump() error {
    dumpFile, err := os.Create("heat-img.dat")
    if err != nil { return err }
    defer func() {
        if err := dumpFile.Close(); err != nil {
            panic(err)
        }
    }()
    // Create a new buffered writer
    w := bufio.NewWriter(dumpFile)
    for i := range s.newMesh.cells {
        for j := range s.newMesh.cells[0] {
            fmt.Fprintf(w, "%f", s.newMesh.cells[i][j])
            if j != len(s.newMesh.cells) - 1 {
                fmt.Fprintf(w, " ")
            }
        }
        fmt.Fprintln(w)
    }
    return w.Flush()
}


func main() {
    var cpuprofile = flag.String("cpuprofile", "", "write CPU profile to file")
    // Parse user input
    flag.Parse()
    // Determine whther or not CPU profiling is on
    if *cpuprofile != "" {
        f, err := os.Create(*cpuprofile)
        if err != nil {
            log.Fatal(err)
        }
        pprof.StartCPUProfile(f)
        defer pprof.StopCPUProfile()
    }
    // Let the games begin
    fmt.Println("o", AppName, AppVerStr)
    sim := NewHeatTxSim(N, N, ThermCond, TMax)
    sim.oldMesh.SetInitConds()
    sim.Run()
    err := sim.Dump()
    if (err != nil) { panic(err) }
}
