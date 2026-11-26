package gnfs

import (
	"bufio"
	"fmt"
	"io"
	"math/big"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

// SaveDir, if set, is a directory where GNFS writes the polynomial and
// relations for each number it factors, as <n>.gnfs.rel.
var SaveDir string

func saveRelations(n *big.Int, poly *Poly, rels []Relation) {
	path := filepath.Join(SaveDir, n.String()+".gnfs.rel")
	f, err := os.Create(path)
	if err == nil {
		err = writeRelations(f, n, poly, rels)
		if cerr := f.Close(); err == nil {
			err = cerr
		}
	}
	if err != nil {
		logf("gnfs: saving relations: %v", err)
		return
	}
	logf("gnfs: saved %d relations to %s", len(rels), path)
}

func writeRelations(w io.Writer, n *big.Int, poly *Poly, rels []Relation) error {
	bw := bufio.NewWriter(w)
	fmt.Fprintln(bw, n)
	fmt.Fprintf(bw, "%d %s %s %s %.17g %.17g\n", poly.Degree, poly.A1, poly.A0, poly.M, poly.Skew, poly.Alpha)
	for i, c := range poly.Coeffs {
		if i > 0 {
			bw.WriteByte(' ')
		}
		bw.WriteString(c.String())
	}
	bw.WriteByte('\n')
	for _, rel := range rels {
		fmt.Fprintf(bw, "%d %d", rel.A, rel.B)
		if rel.IsFree {
			bw.WriteString(" F\n")
			continue
		}
		bw.WriteByte('\n')
	}
	return bw.Flush()
}

func readPolyAndPairs(r io.Reader) (n *big.Int, poly *Poly, pairs [][2]int64, isFree []bool, err error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(nil, 64<<20)
	line := 0
	for sc.Scan() {
		line++
		text := strings.TrimSpace(sc.Text())
		if text == "" {
			continue
		}
		fields := strings.Fields(text)
		switch line {
		case 1:
			var ok bool
			if n, ok = new(big.Int).SetString(text, 10); !ok {
				return nil, nil, nil, nil, fmt.Errorf("line 1: invalid n")
			}
		case 2:
			if len(fields) < 6 {
				return nil, nil, nil, nil, fmt.Errorf("line 2: invalid poly header")
			}
			deg, err := strconv.Atoi(fields[0])
			if err != nil {
				return nil, nil, nil, nil, err
			}
			a1, ok1 := new(big.Int).SetString(fields[1], 10)
			a0, ok0 := new(big.Int).SetString(fields[2], 10)
			m, okm := new(big.Int).SetString(fields[3], 10)
			if !ok1 || !ok0 || !okm {
				return nil, nil, nil, nil, fmt.Errorf("line 2: invalid big.Int")
			}
			skew, _ := strconv.ParseFloat(fields[4], 64)
			alpha, _ := strconv.ParseFloat(fields[5], 64)
			poly = &Poly{Degree: deg, A1: a1, A0: a0, M: m, Skew: skew, Alpha: alpha}
		case 3:
			if len(fields) != poly.Degree+1 {
				return nil, nil, nil, nil, fmt.Errorf("line 3: expected %d coeffs, got %d", poly.Degree+1, len(fields))
			}
			poly.Coeffs = make([]*big.Int, len(fields))
			for i, f := range fields {
				c, ok := new(big.Int).SetString(f, 10)
				if !ok {
					return nil, nil, nil, nil, fmt.Errorf("line 3: invalid coeff %q", f)
				}
				poly.Coeffs[i] = c
			}
		default:
			a, err := strconv.ParseInt(fields[0], 10, 64)
			if err != nil {
				return nil, nil, nil, nil, err
			}
			b, err := strconv.ParseInt(fields[1], 10, 64)
			if err != nil {
				return nil, nil, nil, nil, err
			}
			pairs = append(pairs, [2]int64{a, b})
			isFree = append(isFree, len(fields) > 2 && fields[2] == "F")
		}
	}
	if err := sc.Err(); err != nil {
		return nil, nil, nil, nil, err
	}
	if line < 3 {
		return nil, nil, nil, nil, fmt.Errorf("incomplete relation file")
	}
	return n, poly, pairs, isFree, nil
}
