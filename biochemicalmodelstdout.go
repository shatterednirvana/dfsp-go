package biochemicalstdout

import "fmt"
import "os"

import "./biochemicalmodel"

type BioChemicalSTDOUT struct {
	cnt int
}

func NewBioChemicalSTDOUT() *BioChemicalSTDOUT {
	b := new(BioChemicalSTDOUT)
	b.cnt = 1
	return b
}

func (b BioChemicalSTDOUT) PrintOutput(t float64, m biochemicalmodel.BioChemicalModel, xx [][]int) {
	tspan := m.SampleTimepoints()

	if t == tspan[0] {
		b.cnt = 1
		fmt.Fprint(os.Stdout, "space_arr = -2*pi:", m.VoxelLen(), ":(2*pi-", m.VoxelLen(), ");\n")
		fmt.Fprint(os.Stdout, "sim_output = zeros(", m.NumVoxel(), ",1,", m.NumSpecies(), ");time_arr = zeros(1);\n")
	}

	fmt.Fprint(os.Stdout, "time_arr(", b.cnt, ")=", t, ";\n")

	for s := 0; s < m.NumSpecies(); s++ {
		fmt.Fprint(os.Stdout, "sim_output(:,", b.cnt, ",", s+1, ")=[")
		for i := 0; i < m.NumVoxel(); i++ {
			fmt.Fprint(os.Stdout, " ", xx[i][s])
		}
		fmt.Fprint(os.Stdout, "];\n")
	}
	fmt.Fprint(os.Stderr, t, "\n")

	if t == tspan[len(tspan)-1] {
		for s := 0; s < m.NumSpecies(); s++ {
			fmt.Fprint(os.Stdout, "figure(", s+1, ");surf(time_arr,space_arr,sim_output(:,:,", s+1, "));title('", m.SpeciesName()[s], "');\n")
			fmt.Fprint(os.Stdout, "xlabel('time');ylabel('space'); zlabel('molecules');\n")
		}
	}

	b.cnt++
}
