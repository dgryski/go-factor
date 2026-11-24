package qs

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

// SaveDir, if set, is a directory where the sieve writes the relations for
// each number it factors, as <n>.rel, before solving the matrix.
var SaveDir string

func saveRelations(n *big.Int, factorBase []int64, relations []relation) {
	path := filepath.Join(SaveDir, n.String()+".rel")
	f, err := os.Create(path)
	if err == nil {
		err = writeRelations(f, n, factorBase, relations)
		if cerr := f.Close(); err == nil {
			err = cerr
		}
	}
	if err != nil {
		logf("qs: saving relations: %v", err)
		return
	}
	logf("qs: saved %d relations to %s", len(relations), path)
}

// writeRelations writes n, the factor base and the relations as text. The
// first line is n, the second the factor base, and each further line is a
// relation's X followed by its Factors.
func writeRelations(w io.Writer, n *big.Int, factorBase []int64, relations []relation) error {
	bw := bufio.NewWriter(w)
	fmt.Fprintln(bw, n)
	for i, p := range factorBase {
		if i > 0 {
			bw.WriteByte(' ')
		}
		bw.WriteString(strconv.FormatInt(p, 10))
	}
	bw.WriteByte('\n')
	for _, rel := range relations {
		bw.WriteString(rel.X.String())
		for _, c := range rel.Factors {
			bw.WriteByte(' ')
			bw.WriteString(strconv.Itoa(int(c)))
		}
		bw.WriteByte('\n')
	}
	return bw.Flush()
}

// readRelations reads what writeRelations wrote.
func readRelations(r io.Reader) (n *big.Int, factorBase []int64, relations []relation, err error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(nil, 64<<20)
	line := 0
	for sc.Scan() {
		line++
		fields := strings.Fields(sc.Text())
		switch line {
		case 1:
			var ok bool
			if n, ok = new(big.Int).SetString(sc.Text(), 10); !ok {
				return nil, nil, nil, fmt.Errorf("line 1: invalid n %q", sc.Text())
			}
		case 2:
			for _, f := range fields {
				p, err := strconv.ParseInt(f, 10, 64)
				if err != nil {
					return nil, nil, nil, fmt.Errorf("line 2: %v", err)
				}
				factorBase = append(factorBase, p)
			}
		default:
			if len(fields) == 0 {
				continue
			}
			x, ok := new(big.Int).SetString(fields[0], 10)
			if !ok {
				return nil, nil, nil, fmt.Errorf("line %d: invalid X %q", line, fields[0])
			}
			rel := relation{X: x, Factors: make([]int32, 0, len(fields)-1)}
			for _, f := range fields[1:] {
				c, err := strconv.ParseInt(f, 10, 32)
				if err != nil || c < 0 || int(c) > len(factorBase) {
					return nil, nil, nil, fmt.Errorf("line %d: invalid column %q", line, f)
				}
				rel.Factors = append(rel.Factors, int32(c))
			}
			relations = append(relations, rel)
		}
	}
	if err := sc.Err(); err != nil {
		return nil, nil, nil, err
	}
	if line < 2 {
		return nil, nil, nil, fmt.Errorf("missing n or factor base")
	}
	return n, factorBase, relations, nil
}
