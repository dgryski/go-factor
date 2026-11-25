package factor

// smallPrimes are the primes below 50000, for trial division and p-1.
var smallPrimes = func() []int64 {
	ps := sievePrimes(50000)
	primes := make([]int64, len(ps))
	for i, p := range ps {
		primes[i] = int64(p)
	}
	return primes
}()
