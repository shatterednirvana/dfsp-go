package biochemicalsolver

import "./biochemicalmodeloutput"

type BioChemicalSolver interface {
	Solve(out biochemicalmodeloutput.BioChemicalModelOutput)
}
