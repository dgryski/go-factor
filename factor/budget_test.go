package main

import (
	"math/big"
	"strings"
	"testing"
	"time"
)

func TestBudgetPicksSieve(t *testing.T) {
	n80, _ := new(big.Int).SetString(strings.Repeat("9", 80), 10)
	n105, _ := new(big.Int).SetString(strings.Repeat("9", 105), 10)

	b := &budget{
		start:    time.Now(),
		workers:  15,
		qsRuns:   func(n *big.Int) bool { return n.BitLen() <= qsMaxBits },
		gnfsRuns: func(n *big.Int) bool { return n.BitLen() <= gnfsMaxBits },
	}

	ctx, cancel, name, ok := b.context(n80, shareECM)
	if !ok || name != "qs" {
		t.Fatalf("80 digits got ok=%v name=%q, want qs", ok, name)
	}
	if _, hasDL := ctx.Deadline(); !hasDL {
		t.Fatal("expected deadline for 80 digits")
	}
	cancel()

	ctx, cancel, name, ok = b.context(n105, shareECM)
	if !ok || name != "gnfs" {
		t.Fatalf("105 digits got ok=%v name=%q, want gnfs", ok, name)
	}
	if _, hasDL := ctx.Deadline(); !hasDL {
		t.Fatal("expected deadline for 105 digits")
	}
	cancel()

	nHuge, _ := new(big.Int).SetString(strings.Repeat("9", 150), 10)
	ctx, cancel, name, ok = b.context(nHuge, shareECM)
	if !ok || name != "" {
		t.Fatalf("150 digits got ok=%v name=%q, want no sieve", ok, name)
	}
	if _, hasDL := ctx.Deadline(); hasDL {
		t.Fatal("expected no deadline when neither sieve runs")
	}
	cancel()
}
