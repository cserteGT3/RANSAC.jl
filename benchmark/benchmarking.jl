using Pkg
#Pkg.activate(".")

using Dates
using Logging
using Random

using BenchmarkTools
using CSV
using Tables
using PrettyTables

include(joinpath(dirname(pwd()), "src", "RANSAC.jl"))
using .RANSAC

const table_header = ["date" "commit sha" "minimum time" "median time" "mean time" "maximum time" "allocs" "memory"; "" "" "[s]" "[s]" "[s]" "[s]" "" "[MiB]"]
const nums = :r
const table_align = [:l, :l, nums, nums, nums, nums, nums, nums]

function setupme(iterations)
    vs, ns, norms4Plot, shape_s = examplepc3()
    pcr = PointCloud(vs, ns, 32)
    # plane
    p_ae = (ϵ = 0.3, α=deg2rad(5))
    # cylidner
    cy_ae = (ϵ = 0.3, α=deg2rad(5))
    # sphere
    sp_ae = (ϵ = 0.3, α=deg2rad(5))
    one_ae = AlphSilon(sp_ae, p_ae, cy_ae)
    # number of minimal subsets drawed in one iteration
    tt = 15
    # probability that we found shapes
    ptt = 0.9
    # minimum shape size
    ττ = 900
    # maximum number of iteration
    itermax = iterations
    # size of the minimal set
    draws = 3
    return (pcr, one_ae, tt, ptt, ττ, itermax, draws, 500, true)
end

function makenamedtuple(bmresult)
    comsha = read(`git log -n 1 --pretty=format:"%H"`, String)
    mint = BenchmarkTools.minimum(bmresult).time
    maxt = BenchmarkTools.maximum(bmresult).time
    meant = BenchmarkTools.mean(bmresult).time
    mediant = BenchmarkTools.median(bmresult).time
    alld = allocs(bmresult)
    mmm = memory(bmresult)/1024
    (date=Dates.now(), commitsha=comsha, minimumtime=mint, mediantime=mediant, meantime=meant, maximumtime=maxt, allocated=alld, memory=mmm)
end

function runbenchmark(show = true; savetocsv = false, printlog = false)
    globlog = global_logger()

    prept = setupme(20)
    @info "Precompiling benchmark..."

    printlog && global_logger(NullLogger())
    @benchmark ransac($prept...)
    printlog && global_logger(globlog)

    Random.seed!(1234);
    sharp = setupme(100_000)
    @info "Running real benchmark..."

    printlog && global_logger(NullLogger())
    benched = @benchmark ransac($sharp...)
    printlog && global_logger(globlog)

    @info "Benchmark finished."
    benchedtuple = makenamedtuple(benched)
    bmi = Tables.rowtable([benchedtuple])
    show && printresult(bmi)
    savetocsv && savebenchmark(bmi)
    bmi, benched
end

function nicifyone(onel::Array{T}, divisor) where T<:NamedTuple
    onel = onel[1]
    shas = getproperty(onel, :commitsha)[1:7]
    mint = getproperty(onel, :minimumtime)/divisor
    maxt = getproperty(onel, :maximumtime)/divisor
    meant = getproperty(onel, :meantime)/divisor
    mediant = getproperty(onel, :mediantime)/divisor
    mem = getproperty(onel, :memory)/1024 # KiB to MiB
    [getproperty(onel, :date) shas mint mediant meant maxt getproperty(onel, :allocated) mem]
end

function nicifyone(onel, divisor)
    shas = getproperty(onel, :commitsha)[1][1:7]
    mint = getproperty(onel, :minimumtime)[1]/divisor
    maxt = getproperty(onel, :maximumtime)[1]/divisor
    meant = getproperty(onel, :meantime)[1]/divisor
    mediant = getproperty(onel, :mediantime)[1]/divisor
    mem = getproperty(onel, :memory)[1]/1024 # KiB to MiB
    [getproperty(onel, :date)[1] shas mint mediant meant maxt getproperty(onel, :allocated)[1] mem]
end

function nicify(tb)
    divisor = 1_000_000_000 # nanosec to sec
    if size(tb, 1) < 2
        return nicifyone(tb, divisor)
    end
    shas = [ prop[1:7] for prop in getproperty(tb, :commitsha)]
    mint = [ m/divisor for m in getproperty(tb, :minimumtime)]
    maxt = [ m/divisor for m in getproperty(tb, :maximumtime)]
    meant = [ m/divisor for m in getproperty(tb, :meantime)]
    mediant = [ m/divisor for m in getproperty(tb, :mediantime)]
    mem = [ m/1024 for m in getproperty(tb, :memory)] # KiB to MiB
    [getproperty(tb, :date) shas mint mediant meant maxt getproperty(tb, :allocated) mem]
end

function printresult(tb)
    hmnice = nicify(tb)
    tf = PrettyTableFormat(unicode)
    pretty_table(hmnice, table_header, tf, alignment = table_align, formatter=ft_round(3, [3,4,5,6,8]))
end

function saveresult(tb)
    hmnice = nicify(tb)
    tf = PrettyTableFormat(markdown)
    fname = "benchmark_results.md"
    open(fname, "w") do io
        pretty_table(io, hmnice, table_header, tf, alignment = table_align, formatter=ft_round(3, [3,4,5,6,8]))
    end
    @info "Saved file."
end

function savebenchmark(bm)
    fname = "benchmark_results.csv"
    if !isfile(fname)
        CSV.write(fname, bm )
    else
        CSV.write(fname, bm, append = true)
    end
    @info "Saved file."
end

function loadbenchmarks()
    fname = "benchmark_results.csv"
    return CSV.read(fname)
end

function info()
    @info "Runnnig one benchmark takes around 2-3 minutes."
    @info "`bmark, benched = runbenchmark(true)` to run the benchmark (and also display it)."
    @info "`savebenchmark(bmark)` to append the last benchmark to the CSV file."
    @info "`loadbenchmarks()` to load the saved benchmarks."
    @info "`printresult(bmark)` to show the result of one or more benchmarks."
    @info "`saveresult(bmark)` to save the prettyprint to markdown."
end

info()
