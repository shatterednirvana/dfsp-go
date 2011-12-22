package dfspsolver

import "fmt"
import "math"
import "os"
import "rand"

import "./biochemicalmodeloutput"
import "./gproteincycle"
import "./matrix"

type DFSPSolver struct {
	TauD float64
	DiffuseErrTol float64
	MaxJumps int

	SystemState [][]int
	ARxn [][]float64

	CacheInitialized bool
	Cache []float64
	Model gproteincycle.GProteinCycle
	//model biochemicalmodel.BioChemicalModel

	SampleTimepoints int
	NumRxn int
	NumDiff int	
}

func NewDFSPSolver(model gproteincycle.GProteinCycle) *DFSPSolver {
	d := new(DFSPSolver)
	d.MaxJumps = 5
	isSetTau := false
	d.CacheInitialized = false
	d.Model = model

	if !isSetTau {
		d.TauD = 0.1
	}

	d.ARxn = make([][]float64, d.Model.NumVoxel())
	for i := 0; i < d.Model.NumVoxel(); i++ {
		d.ARxn[i] = make([]float64, d.Model.NumRxn())
	}

	d.NumRxn = 0
	d.NumDiff = 0
	return d
}

func (d DFSPSolver) DumpA() {
	for j := 0; j < d.Model.NumRxn(); j++ {
		fmt.Fprint(os.Stdout, "\t", j, ":")
		for i := 0; i < d.Model.NumVoxel(); i++ {
			fmt.Fprint(os.Stdout, d.ARxn[i][j], " ")
		}
	}
}

func (d *DFSPSolver) Solve(output biochemicalmodeloutput.BioChemicalModelOutput) {
	t := 0.0
	tFinal := d.Model.SampleTimepoints()[len(d.Model.SampleTimepoints())-1]
	d.SystemState = d.Model.InitializeState()
	tNdx := 0
	tNextD := 0.0
	a0 := d.InitializeRxnPropensity()
	fmt.Fprint(os.Stderr, "a0=", a0)
	d.DumpA()
	tNextR := 0.0

	if a0 > 0 {
		tNextR = t - math.Log(rand.Float64()) / a0
	} else {
		tNextR = tFinal + 1
	}

	fmt.Fprint(os.Stderr, "t_next_R", tNextR)
	tNextO := d.Model.SampleTimepoints()[tNdx]
	d.SampleTimepoints = 0
	d.NumRxn = 0
	d.NumDiff  = 0

	//for t < tFinal {
	for t < 1 {
		if tNextO <= tNextD && tNextO <= tNextR {
			t = d.Model.SampleTimepoints()[tNdx]
			output.PrintOutput(t, d.Model, d.SystemState)
			tNdx++
			if tNdx >= len(d.Model.SampleTimepoints()) {
				break
			} else {
				tNextO = d.Model.SampleTimepoints()[tNdx]
			}
			d.SampleTimepoints++
		} else if tNextD <= tNextO && tNextD <= tNextR {
			t = tNextD
			d.TakeDiffusionStep()
			tNextD = t + d.TauD
			a0 = d.InitializeRxnPropensity()
			if a0 > 0 {
				tNextR = t - math.Log(rand.Float64()) / a0
			} else {
				tNextR = tFinal + 1
			}
			d.NumDiff++
		} else if tNextR <= tNextO && tNextR <= tNextD {
			t = tNextR
			a0 = d.TakeRxnStep(a0)
			if a0 > 0 {
				tNextR = t - math.Log(rand.Float64()) / a0
			} else {
				tNextR = tFinal + 1
			}
			d.NumRxn++
		} else {
			fmt.Fprint(os.Stderr, "DFSP: should not get here")
			os.Exit(1)
		}
	}
}

func (d DFSPSolver) PrintSimStats() {
	fmt.Fprint(os.Stderr, "# of s timepoints = ", d.SampleTimepoints, "\n")
	fmt.Fprint(os.Stderr, "# of reactions = ", d.NumRxn, "\n")
	fmt.Fprint(os.Stderr, "# of diffusions = ", d.NumDiff, "\n")
}

func (d *DFSPSolver) TakeDiffusionStep() {
	nextX := make([][]int, d.Model.NumVoxel())
	for i := 0; i < d.Model.NumVoxel(); i++ {
		nextX[i] = make([]int, d.Model.NumSpecies())
	}

	for s := 0; s < d.Model.NumSpecies(); s++ {
		if d.Model.DiffusionCoefficients()[s] > 0 {
			for i := 0; i < d.Model.NumVoxel(); i++ {
				subX := d.TakeStepFsp(d.SystemState[i][s], d.TauD, d.Model.DiffusionCoefficients()[s], d.Model.VoxelLen())
				fst := i - int(math.Floor(float64(len(subX))/float64(2.0)))
				for j := 0; j < len(subX); j++ {
					k := fst+j
					if k >= d.Model.NumVoxel() {
						k -= d.Model.NumVoxel()
					}

					if k < 0 {
						k += d.Model.NumVoxel()
					}

					nextX[k][s] += subX[j]
				}
			}
		} else {			
			for i := 0; i < d.Model.NumVoxel(); i++ {
				nextX[i][s] = d.SystemState[i][s]
			}
		}
	}
	d.SystemState = nextX
}

func (d DFSPSolver) TakeStepFsp(n int, tau float64, D float64, l float64) []int {
	if !d.CacheInitialized {
		a := d.GenerateAMatrix(1)
		tauParam := tau * D / (l*l)
		w := d.CalculateProbabilities(a, tauParam)
		d.CacheInitialized = true
		d.Cache = w
	}

	out := make([]int, 2*d.MaxJumps+1)
	for i := 0; i < n; i++ {
		rnd := rand.Float64()
		sum := 0.0
		ndx := 0
		for k := 0; k < len(d.Cache); k++ {
			sum += d.Cache[k]
			if sum > rnd {
				ndx = k
				break
			}
		}
		finalPos := d.GenerateStates(1, ndx)
		for j := 0; j < 2*d.MaxJumps+1; j++ {
			out[j] += finalPos[j]
		}
	}
	return out
}

func (d DFSPSolver) CalculateProbabilities(a [][]int, tau float64) []float64 {
	m := matrix.NewMatrixFromIntData(a)
	m = m.TimesScalar(tau)
	m = m.Exponential()
	wT := m.GetColumn(0)
	w := make([]float64, len(wT)-1)
	for i := 0; i < len(w); i++ {
		w[i] = wT[i] / (1.0 - wT[len(wT)-1])
	}
	return w
}

func (d DFSPSolver) GenerateAMatrix(nI int) [][]int {
	if nI == 1 {
		return [][]int{
			{-2,1,1,0,0,0,0,0,0,0,0,0},
			{1,-2,0,1,0,0,0,0,0,0,0,0},
			{1,0,-2,0,1,0,0,0,0,0,0,0},
			{0,1,0,-2,0,1,0,0,0,0,0,0},
			{0,0,1,0,-2,0,1,0,0,0,0,0},
			{0,0,0,1,0,-2,0,1,0,0,0,0},
			{0,0,0,0,1,0,-2,0,1,0,0,0},
			{0,0,0,0,0,1,0,-2,0,1,0,0},
			{0,0,0,0,0,0,1,0,-2,0,1,0},
			{0,0,0,0,0,0,0,1,0,-2,0,1},
			{0,0,0,0,0,0,0,0,1,0,-2,1},
			{0,0,0,0,0,0,0,0,0,1,1,0}}
	}

	fmt.Fprint(os.Stderr, "DFSP: not implemented yet\n")
	os.Exit(1)

	return [][]int{
		{-2,1,1,0,0,0,0,0,0,0,0,0},
		{1,-2,0,1,0,0,0,0,0,0,0,0},
		{1,0,-2,0,1,0,0,0,0,0,0,0},
		{0,1,0,-2,0,1,0,0,0,0,0,0},
		{0,0,1,0,-2,0,1,0,0,0,0,0},
		{0,0,0,1,0,-2,0,1,0,0,0,0},
		{0,0,0,0,1,0,-2,0,1,0,0,0},
		{0,0,0,0,0,1,0,-2,0,1,0,0},
		{0,0,0,0,0,0,1,0,-2,0,1,0},
		{0,0,0,0,0,0,0,1,0,-2,0,1},
		{0,0,0,0,0,0,0,0,1,0,-2,1},
		{0,0,0,0,0,0,0,0,0,1,1,0}}
}

func (d DFSPSolver) GenerateStates(nI int, ndx int) [11]int {
	if nI == 1 {
		if ndx == 0 {
			return [...]int{0,0,0,0,0,1,0,0,0,0,0}
		} else if ndx == 1 {
			return [...]int{0,0,0,0,0,0,1,0,0,0,0}
		} else if ndx == 2 {
			return [...]int{0,0,0,0,1,0,0,0,0,0,0}
		} else if ndx == 3 {
			return [...]int{0,0,0,0,0,0,0,1,0,0,0}
		} else if ndx == 4 {
			return [...]int{0,0,0,1,0,0,0,0,0,0,0}
		} else if ndx == 5 {
			return [...]int{0,0,0,0,0,0,0,0,1,0,0}
		} else if ndx == 6 {
			return [...]int{0,0,1,0,0,0,0,0,0,0,0}
		} else if ndx == 7 {
			return [...]int{0,0,0,0,0,0,0,0,0,1,0}
		} else if ndx == 8 {
			return [...]int{0,1,0,0,0,0,0,0,0,0,0}
		} else if ndx == 9 {
			return [...]int{0,0,0,0,0,0,0,0,0,0,1}
		} else if ndx == 10 {
			return [...]int{1,0,0,0,0,0,0,0,0,0,0}
		}
	} else {
		fmt.Fprint(os.Stderr, "DFSP: not implemented yet")
		os.Exit(1)
	}

	fmt.Fprint(os.Stderr, "should not get here")
	os.Exit(1)

	return [...]int{1,0,0,0,0,0,0,0,0,0,0}
}

func (d *DFSPSolver) TakeRxnStep(a0 float64) float64 {
	c := rand.Float64() * a0
	tmp := 0.0
	i := 0
	r := 0
	for i := 0; i < d.Model.NumVoxel(); i++ {
		for r := 0; r < d.Model.NumRxn(); r++ {
			tmp += d.ARxn[i][r]
			if tmp > c {
				break
			}
		}
		if tmp > c {
			break
		}
	}
	d.Model.DoReaction(i, r, d.SystemState)
	for r := 0; r < d.Model.NumRxn(); r++ {
		a0 -= d.ARxn[i][r]
		d.ARxn[i][r] = d.Model.UpdateRxn(i, r, d.SystemState)
		a0 += d.ARxn[i][r]
	}
	return a0
}

func (d *DFSPSolver) InitializeRxnPropensity() float64 {
	a0 := 0.0
	for i := 0; i < d.Model.NumVoxel(); i++ {
		for r := 0; r < d.Model.NumRxn(); r++ {
			d.ARxn[i][r] = d.Model.UpdateRxn(i, r, d.SystemState)
			a0 += d.ARxn[i][r]
		}
	}
	return a0
}
