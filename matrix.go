package matrix

import "fmt"
import "math"
import "os"
import "rand"

type Matrix struct {
	M int
	N int
	Data [][]float64
}

func NewMatrixFromMandN(M int, N int) *Matrix {
	m := new(Matrix)
	m.M = M
	m.N = N
	// TODO - does this initialize Data to zero?
	m.Data = make([][]float64, M)
	for i := 0; i < M; i++ {
		m.Data[i] = make([]float64, N)
	}
	return m
}

func NewMatrixFromIntData(data [][]int) *Matrix {
	M := len(data)
	N := len(data[0])
	
	m := new(Matrix)
	m.M = M
	m.N = N
	m.Data = make([][]float64, M)
	for i := 0; i < M; i++ {
		m.Data[i] = make([]float64, N)
		for j := 0; j < N; j++ {
			m.Data[i][j] = float64(data[i][j])
		}
	}
	return m
}

func NewMatrixFromFloatData(data [][]float64) *Matrix {
	M := len(data)
	N := len(data[0])
	
	m := new(Matrix)
	m.M = M
	m.N = N
	for i := 0; i < M; i++ {
		for j := 0; j < N; j++ {
			m.Data[i][j] = data[i][j]
		}
	}
	return m
}

func NewMatrixFromMatrix(A Matrix) *Matrix {
	m := new(Matrix)
	m.M = A.M
	m.N = A.N
	m.Data = A.Data
	return m
}

func (m Matrix) Size() [2]int {
	return [...]int{m.M, m.N}
}

func Random(m int, n int) *Matrix {
	a := NewMatrixFromMandN(m, n)
	for i := 0; i < m; i++ {
		for j := 0; j < n; j++ {
			a.Data[i][j] = rand.Float64()
		}
	}
	return a
}

func Identity(n int) *Matrix {
	identity := NewMatrixFromMandN(n, n)
	for i := 0; i < n; i++ {
		identity.Data[i][i] = 1
	}
	return identity
}

func (m *Matrix) Swap(i int, j int) {
	temp := (*m).Data[i]
	(*m).Data[i] = (*m).Data[j]
	(*m).Data[j] = temp
}

func (m Matrix) Transpose() *Matrix {
	a := NewMatrixFromMandN(m.N, m.M)
	for i := 0; i < m.M; i++ {
		for j := 0; j < m.N; j++ {
			a.Data[j][i] = m.Data[i][j]
		}
	}
	return a
}

func (m Matrix) Plus(b Matrix) *Matrix {
	a := m
	if b.M != a.M || b.N != a.N {
		fmt.Fprintf(os.Stderr, "Illegal matrix dimensions.")
		os.Exit(1)
	}

	c := NewMatrixFromMandN(a.M, a.N)
	for i := 0; i < a.M; i++ {
		for j := 0; j < a.N; j++ {
			c.Data[i][j] = a.Data[i][j] + b.Data[i][j]
		}
	}
	return c
}

func (m Matrix) Minus(b Matrix) *Matrix {
	a := m
	if b.M != a.M || b.N != a.N {
		fmt.Fprintf(os.Stderr, "Illegal matrix dimensions.")
		os.Exit(1)
	}

	c := NewMatrixFromMandN(a.M, a.N)
	for i := 0; i < a.M; i++ {
		for j := 0; j < a.N; j++ {
			c.Data[i][j] = a.Data[i][j] - b.Data[i][j]
		}
	}
	return c
}

func (m Matrix) Eq(b Matrix) bool {
	a := m
	if b.M != a.M || b.N != a.N {
		fmt.Fprintf(os.Stderr, "Illegal matrix dimensions.")
		os.Exit(1)
	}

	for i := 0; i < a.M; i++ {
		for j := 0; j < a.N; j++ {
			if a.Data[i][j] != b.Data[i][j] {
				return false
			}
		}
	}
	return true
}

func (m Matrix) Times(b Matrix) *Matrix {
	a := m
	if a.N != b.M {
		fmt.Fprintf(os.Stderr, "Illegal matrix dimensions.")
		os.Exit(1)
	}

	c := NewMatrixFromMandN(a.M, b.N)
	for i := 0; i < c.M; i++ {
		for j := 0; j < c.N; j++ {
			sum := float64(0)
			aRow := a.Data[i][0:a.N]
			bRow := b.Data[i][0:a.N]
			for k := 0; k < a.N; k++ {
				sum += (aRow[k] * bRow[k])
			}
			c.Data[i][j] = sum
		}
	}
	return c
}

func (m Matrix) TimesScalar(x float64) *Matrix {
	a := m
	c := NewMatrixFromMandN(a.M, a.N)
	for i := 0; i < a.M; i++ {
		aRow := a.Data[i][0:a.N]
		for j := 0; j < a.N; j++ {
			c.Data[i][j] = aRow[j] * x
		}
	}
	return c
}

// TODO - this could use a slice perhaps?
func (m Matrix) GetColumn(j int) []float64 {
	r := make([]float64, m.M)
	for i := 0; i < m.M; i++ {
		r[i] = m.Data[i][j]
	}
	return r
}

func (m Matrix) Solve(rhs Matrix) *Matrix {
	if m.M != m.N || rhs.M != m.N || rhs.N != 1 {
		fmt.Fprintf(os.Stderr, "Illegal matrix dimensions.")
		os.Exit(1)
	}

	a := NewMatrixFromMatrix(m)
	b := NewMatrixFromMatrix(rhs)

	for i := 0; i < m.N; i++ {
		max := i
		for j := i+1; j < m.N; j++ {
			if math.Fabs(a.Data[j][i]) > math.Fabs(a.Data[max][i]) {
				max = j
			}
		}

		a.Swap(i, max)
		b.Swap(i, max)

		if a.Data[i][i] == 0.0 {
			fmt.Fprintf(os.Stderr, "Matrix is singular.")
			os.Exit(1)
		}

		for j := i+1; j < m.N; j++ {
			b.Data[j][0] -= b.Data[i][0] * a.Data[j][i] / a.Data[i][i]
		}

		for j := i+1; j < m.N; j++ {
			temp := a.Data[j][i] / a.Data[i][i]
			for k := i+1; k < m.N; k++ {
				a.Data[j][k] -= a.Data[i][k] * temp
			}
			a.Data[j][i] = 0.0
		}
	}

	x := NewMatrixFromMandN(m.N, 1)
	for j := m.N-1; j >= 0; j-- {
		t := 0.0
		for k := j+1; k < m.N; k++ {
			t += a.Data[j][k] * x.Data[k][0]
		}
		x.Data[j][0] = (b.Data[j][0] - t) / a.Data[j][j]
	}
	return x
}

func (m Matrix) Exponential() *Matrix {
	a := m

	if a.M != a.N {
		fmt.Fprintf(os.Stderr, "Illegal matrix dimensions, must be square.")
		os.Exit(1)
	}

	mLn2 := 0.693147180559945
	sugg := a.ExponentialObtainSuggestion()
	divisor := math.Exp(mLn2 * float64(sugg[1]))
	reducedA := a.TimesScalar(1/divisor)
	eA := reducedA.MatrixExpSeries(sugg[0])
	for i := 0; i < sugg[1]; i++ {
		eA = (*eA).Times(*eA)
	}
	return eA
}

func (m Matrix) MatrixExpSeries(numTerms int) *Matrix {
	b := m
	eB := b.TimesScalar(1.0/float64(numTerms)).Plus(*Identity(b.M))
	for count := numTerms-1; count >= 1; count-- { // TODO - original did --count, is that ok?
		eB = eB.Times(b).TimesScalar(1.0/float64(count)).Plus(*Identity(b.M))
	}
	return eB
}

func (m Matrix) ExponentialObtainSuggestion() [2]int {
	normA := m.SupNorm()
	if normA < 0.01 {
		return [...]int{5,1}
	} else if normA < 0.1 {
		return [...]int{5, 4}
	} else if normA < 1.0 {
		return [...]int{7, 5}
	} else if normA < 10.0 {
		return [...]int{9, 7}
	} else if normA < 100.0 {
		return [...]int{10, 10}
	} else if normA < 1000.0 {
		return [...]int{8, 14}
	} else {
		fmt.Fprintf(os.Stderr, "Matrix exponential for matrix with norm >= 1000.")
		os.Exit(1)
	}

	return [...]int{0,0}
}

func (m Matrix) SupNorm() float64 {
	a := m
	max := 0.0
	for i := 0; i < a.M; i++ {
		for j := 0; j < a.N; j++ {
			newMax := math.Fabs(a.Data[i][j])
			if newMax > max {
				max = newMax
			}
		}
	}
	return max
}

func (m Matrix) Show() {
	for i := 0; i < m.M; i++ {
		for j := 0; j < m.N; j++ {
			fmt.Fprintf(os.Stdout, "%9.4f ", m.Data[i][j])
		}
		fmt.Fprint(os.Stdout, "\n")
	}
}
