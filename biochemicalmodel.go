package biochemicalmodel

type BioChemicalModel interface {
	NumRxn() int
	NumSpecies() int
	NumVoxel() int
	VoxelLen() float64

	SpeciesName() [7]string

	ConnectionMatrix() [][]float64
	DiffusionCoefficients() [7]float64
	StoicMatrix() [8][7]int
	RxnRates() [8]float64
	RxnReagents() [8][7]int

	SampleTimepoints() [51]float64

	InitializeState() [][]int

	//GetArguments() []string
}
