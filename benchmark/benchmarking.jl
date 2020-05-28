module BenchmarkRANSAC

using Dates
using Logging
using Random

using BenchmarkTools
using Tables
using PrettyTables

using RANSAC

export runbenchmark, savebenchmark
export loadbenchmarks
export printresult, saveresult
export info
export printload

const table_header = ["date" "commit sha" "minimum time" "median time" "mean time" "maximum time" "allocs" "memory" "system"; "" "" "[s]" "[s]" "[s]" "[s]" "" "[MiB]" ""]
const md_table_header = ["date" "commit sha" "minimum time [s]" "median time [s]" "mean time [s]" "maximum time [s]" "allocs" "memory [MiB]" "system"]
const nums = :r
const table_align = [:l, :l, nums, nums, nums, nums, nums, nums, :l]

function setupme(iterations)
    vs, ns, norms4Plot, shape_s = examplepc3()
    pcr = RANSACCloud(vs, ns, 32)
    #p = RANSACParameters{Float64}()
    # plane
    # p = RANSACParameters(p, ϵ_plane=0.3, α_plane=deg2rad(5))
    p_plane = (ϵ=0.3, α=deg2rad(5),)
    # cylidner
    # p = RANSACParameters(p, ϵ_cylinder=0.3, α_cylinder=deg2rad(5))
    p_cyl = (ϵ=0.3, α=deg2rad(5),)
    # sphere
    # p = RANSACParameters(p, ϵ_sphere=0.3, α_sphere=deg2rad(5))
    p_sp = (ϵ=0.3, α=deg2rad(5),)
    # number of minimal subsets drawed in one iteration
	# probability that we found shapes
	# minimum shape size
	# maximum number of iteration
    p_it = (minsubsetN=15, prob_det=0.9, τ=900, itermax=iterations, drawN=3)
    p = ransacparameters(;plane=p_plane, cylinder=p_cyl, sphere=p_sp, iteration=p_it)
    return (pcr, p, true)
end

function makenamedtuple(bmresult)
    comsha = read(`git log -n 1 --pretty=format:"%H"`, String)
    pc = getpc()
    mint = BenchmarkTools.minimum(bmresult).time
    maxt = BenchmarkTools.maximum(bmresult).time
    meant = BenchmarkTools.mean(bmresult).time
    mediant = BenchmarkTools.median(bmresult).time
    alld = allocs(bmresult)
    mmm = memory(bmresult)/1024
    (date=Dates.now(), commitsha=comsha, minimumtime=mint, mediantime=mediant, meantime=meant, maximumtime=maxt, allocated=alld, memory=mmm, system = pc)
end

function getpc()
    if Sys.iswindows()
        username = ENV["UserName"]
        if username == "cstamas"
            return "WorkLaptop"
        elseif username == "Ipse"
            return "HomeLaptop"
		elseif username == "Laci"
            return "HomePC"
        else
            return "unknownWin"
        end
    else
        return "notWin"
    end
end

"""
    runbenchmark(show = true; savetocsv = false, showdebug = false)

Run benchmark.
"""
function runbenchmark(show = true; savetocsv = false, showdebug = false)
    glb = global_logger()
    showdebug ? global_logger(nosource_debuglogger()) : global_logger(nosource_infologger())

    Random.seed!(1234);
    sharp = setupme(100_000)
    @info "Running real benchmark..."
    benched = @benchmark ransac($sharp..., reset_rand=true)
    @info "Benchmark finished."

    benchedtuple = makenamedtuple(benched)
    bmi = Tables.rowtable([benchedtuple])
    show && display(benched)
    savetocsv && savebenchmark(bmi)
    global_logger(glb)
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
    sys = getproperty(onel, :system)
    [getproperty(onel, :date) shas mint mediant meant maxt getproperty(onel, :allocated) mem sys]
end

function nicifyone(onel, divisor)
    shas = getproperty(onel, :commitsha)[1][1:7]
    mint = getproperty(onel, :minimumtime)[1]/divisor
    maxt = getproperty(onel, :maximumtime)[1]/divisor
    meant = getproperty(onel, :meantime)[1]/divisor
    mediant = getproperty(onel, :mediantime)[1]/divisor
    mem = getproperty(onel, :memory)[1]/1024 # KiB to MiB
    sys = getproperty(onel, :system)
    [getproperty(onel, :date)[1] shas mint mediant meant maxt getproperty(onel, :allocated)[1] mem sys]
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
    sys = getproperty(tb, :system)
    [getproperty(tb, :date) shas mint mediant meant maxt getproperty(tb, :allocated) mem sys]
end

function printresult(tb)
    hmnice = nicify(tb)
    tf = PrettyTableFormat(unicode)
    pretty_table(hmnice, table_header, tf, alignment = table_align, formatter=ft_round(3, [3,4,5,6,8]))
end

function info()
    @info "Runnnig one benchmark takes around 5 minutes."
    @info "`bmark, benched = runbenchmark();` to run the benchmark (and also display it)."
    @info "`printresult(bmark)` to show the result of one or more benchmarks."
end

end # module

using .BenchmarkRANSAC
info()
