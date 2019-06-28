# Benchmarks

Reason of this folder is to be able to benchmark the current state of the code for comaparison.
This means that the I have to save and load the results.
Also plotting would be nice, but first a table is fine.
Benchmarks must be identifiable, I use commit hash and date for this purpose.

Benchmark must be run from a different session and there should be a switch for saving/not saving the results.

## Saving

The saved `.csv` file should look like (every time in nanoseconds, memory in KiB):
```
date, commit sha, minimumtime, mediantime, meantime, maximumtime, allocs, memory
```

Commit hash can be retreived with (first seven character is the short hash on GitHub):
```
read(`git log -n 1 --pretty=format:"%H"`, String)
```

## Interface

~~It uses the repo's root project file (if present).~~
Project file will be used if [27418](https://github.com/JuliaLang/julia/issues/27418) will be solved.

The following packages are needed to run and print the benchmarks:
* [BenchmarkTools](https://github.com/JuliaCI/BenchmarkTools.jl)
* [CSV](https://github.com/JuliaData/CSV.jl)
* [PrettyTables](https://github.com/ronisbr/PrettyTables.jl)
* [Tables](https://github.com/JuliaData/Tables.jl)

```
> cd root/to/repo/benchmarks

> runbenchmarks.bat
[ Info: Runnnig one benchmark takes around 2-3 minutes.
[ Info: `bmark, benched = runbenchmark(true)` to run the benchmark (and also display it).
[ Info: `savebenchmark(bmark)` to append the last benchmark to the CSV file.
[ Info: `loadbenchmarks()` to load the saved benchmarks.
[ Info: `printresult(bmark)` to show the result of one or more benchmarks.
[ Info: `saveresult(bmark)` to save the prettyprint to markdown.
```
