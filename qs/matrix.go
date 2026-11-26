package qs

import "github.com/dgryski/go-factor/internal/linalg"

const keepSurplus = linalg.KeepSurplus

func findDependencies(relations []relation, numCols int) [][]uint64 {
	rows := make([][]int32, len(relations))
	for i, rel := range relations {
		rows[i] = linalg.OddColumns(rel.Factors)
	}
	return linalg.FindDependencies(rows, numCols, func(f string, a ...interface{}) {
		logf("qs: "+f, a...)
	})
}
