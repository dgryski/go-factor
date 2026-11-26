package main

import (
	"math/big"
	"strings"
	"testing"
	"time"
)

func TestBudgetPicksSieve(t *testing.T) {
	n70, _ := new(big.Int).SetString(strings.Repeat("9", 70), 10)
	n85, _ := new(big.Int).SetString(strings.Repeat("9", 85), 10)

	b := &budget{
		start:    time.Now(),
		workers:  15,
		qsRuns:   func(n *big.Int) bool { return n.BitLen() <= 266 },
		gnfsRuns: func(n *big.Int) bool { return n.BitLen() <= 450 },
	}

	ctx, cancel, name, ok := b.context(n70, shareECM)
	if !ok || name != "qs" {
		t.Fatalf("70 digits got ok=%v name=%q, want qs", ok, name)
	}
	if _, hasDL := ctx.Deadline(); !hasDL {
		t.Fatal("expected deadline for 70 digits")
	}
	cancel()

	ctx, cancel, name, ok = b.context(n85, shareECM)
	if !ok || name != "gnfs" {
		t.Fatalf("85 digits got ok=%v name=%q, want gnfs", ok, name)
	}
	if _, hasDL := ctx.Deadline(); !hasDL {
		t.Fatal("expected deadline for 85 digits")
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
