package gproteincycle

import "fmt"
import "math"
import "os"

type GProteinCycle struct {

}

func (g GProteinCycle) NumRxn() int {
	return 8
}

func (g GProteinCycle) NumSpecies() int {
	return 7
}

func (g GProteinCycle) NumVoxel() int {
	return 200
}

func (g GProteinCycle) VoxelLen() float64 {
	return 0.0628318
}

func (g GProteinCycle) SpeciesName() [7]string {
	return [...]string{"L","R","RL","G","Ga","Gbg","Gd"}
}

func (g GProteinCycle) ConnectionMatrix() [][]float64 {
	arr := make([][]float64, 1)
	arr[0] = make([]float64, 1)
	arr[0][0] = 1
	return arr
}

func (g GProteinCycle) DiffusionCoefficients() [7]float64 {
	return [...]float64{0,0.001,0.001,0.001,0.001,0.001,0.001}
}

func (g GProteinCycle) StoicMatrix() [8][7]int {
	return [...][...]int{{0,1,0,0,0,0,0},
		{0,-1,0,0,0,0,0},
		{0,-1,1,0,0,0,0},
		{0,1,-1,0,0,0,0},
		{0,0,-1,0,0,0,0},
		{0,0,0,-1,1,1,0},
		{0,0,0,0,-1,0,1},
		{0,0,0,1,0,-1,-1}}
}

func (g GProteinCycle) RxnRates() [8]float64 {
	return [...]float64{3.97899,4e-04,1e-7,1e-02,4e-04,1.00528e-6,0.1,1.00528}
}

func (g GProteinCycle) RxnReagents() [8][7]int {
	return [...][...]int{{0,0,0,0,0,0,0},
		{0,1,0,0,0,0,0},
		{1,1,0,0,0,0,0},
		{0,0,1,0,0,0,0},
		{0,0,1,0,0,0,0},
		{0,0,1,1,0,0,0},
		{0,0,0,0,1,0,0},
		{0,0,0,0,0,1,1}}
}

func (g GProteinCycle) SampleTimepoints() [51]float64 {
	return [...]float64{0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100}
}

func (g GProteinCycle) InitializeState() [][]int {
	rt := 10000
	gt := 10000
	lMin := 0
	lMax := 400

	x := make([][]int, g.NumVoxel())
	for i := 0; i < g.NumVoxel(); i++ {
		x[i] = make([]int, g.NumSpecies())
	}

	for i := 0; i < g.NumVoxel(); i++ {
		for j := 0; j < g.NumSpecies(); j++ {
			if j == 0 {
				x[i][j] = int(math.Floor(float64(lMax-lMin)*0.5*float64(1+math.Cos(0.5*(float64(i)*float64(g.VoxelLen())-2*3.14159)))+float64(lMin)))
			} else if j == 1 {
				x[i][j] = int(math.Floor(float64(rt)/float64(g.NumVoxel())))
			} else if j == 3 {
				x[i][j] = int(math.Floor(float64(gt)/float64(g.NumVoxel())))
			} else {
				x[i][j] = 0
			}
		}
	}

	return x
}

func NewGProteinCycle() *GProteinCycle {
	return new(GProteinCycle)
}

func (g GProteinCycle) UpdateRxn(voxelIndex int, reactionIndex int, systemState [][]int) float64 {
	vol := math.Pow(g.VoxelLen(), 3)
	a := g.RxnRates()[reactionIndex]
	rxnOrder := 0

	for k := 0; k < g.NumSpecies(); k++ {
		if g.RxnReagents()[reactionIndex][k] == 1 {
			rxnOrder += 1
			a *= float64(systemState[voxelIndex][k])
		} else if g.RxnReagents()[reactionIndex][k] == 2 {
			if systemState[voxelIndex][k] >= 2 {
				rxnOrder += 2
				a *= (float64(systemState[voxelIndex][k]) * float64(systemState[voxelIndex][k]-1)) / 2.0
			} else {
				a = 0
				rxnOrder = 0
			}
		}
	}

	if rxnOrder == 0 {
		return a * vol
	} else if rxnOrder == 1 {
		return a
	} else if rxnOrder == 2 {
		return a / vol
	} else {
		fmt.Fprint(os.Stderr, "tri-molecular or above reactions not supported")
		os.Exit(1)
	}

	return 0
}

func (g GProteinCycle) DoReaction(voxelIndex int, reactionIndex int, systemState [][]int) {
	for j := 0; j < g.NumSpecies(); j++ {
		systemState[voxelIndex][j] += g.StoicMatrix()[reactionIndex][j]
	}
}

func (g GProteinCycle) UpdateDiffusion(voxelIndex int, speciesIndex int, systemState [][]int) float64 {
	return float64(g.DiffusionCoefficients()[speciesIndex])*float64(systemState[voxelIndex][speciesIndex]) / float64(g.VoxelLen()*g.VoxelLen())
}

// TODO - I think the array is passed by value and not reference - fix that
func (g GProteinCycle) DoDiffusion(srcVoxelIndex int, destVoxelIndex int, speciesIndex int, systemState [][]int) {
	systemState[srcVoxelIndex][speciesIndex]--
	systemState[destVoxelIndex][speciesIndex]++
}
