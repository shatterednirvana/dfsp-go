package main

import "flag"
import "log"
import "os"
import "runtime/pprof"

import "./biochemicalmodelstdout"
import "./dfspsolver"
import "./gproteincycle"

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")

func main() {
	flag.Parse()
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

        g := *gproteincycle.NewGProteinCycle()
        s := dfspsolver.NewDFSPSolver(g)
        s.Solve(biochemicalstdout.NewBioChemicalSTDOUT())
        s.PrintSimStats()
}
