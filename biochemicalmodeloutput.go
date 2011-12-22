package biochemicalmodeloutput

import "./biochemicalmodel"

type BioChemicalModelOutput interface {
	PrintOutput(t float64, m biochemicalmodel.BioChemicalModel, system_state [][]int)
}
