package factor

import "time"

// Logf, if set, receives progress messages from the factoring methods.
var Logf func(format string, args ...interface{})

func logf(format string, args ...interface{}) {
	if Logf != nil {
		Logf(format, args...)
	}
}

func since(t time.Time) time.Duration {
	return time.Since(t).Round(time.Millisecond)
}
