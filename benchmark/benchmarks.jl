using BenchmarkTools
using Random
using RANSAC

const SUITE = BenchmarkGroup()

SUITE["RANSAC"] = BenchmarkGroup()

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

Random.seed!(1234);
sharp = setupme(100_000)
SUITE["RANSAC"]["smallideal"] = @benchmarkable ransac($sharp..., reset_rand=true)