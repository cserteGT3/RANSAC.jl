# Benchmarks

Reason of this folder is to be able to benchmark the current state of the code for comparison.
The `benchmarks.jl` file contains code for [`PkgBenchmark.jl`](https://github.com/JuliaCI/PkgBenchmark.jl).
Old code is in `benchmarking.jl` but it partially stays, because it makes it easy to run the full algorithm without any setup.
Another reason is that with the `showdebug` keyword logging can be adjusted, which I wasn't able with `PkgBenchmark`.
