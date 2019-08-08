using BenchmarkTools
using Random
using RANSAC

const SUITE = BenchmarkGroup()

SUITE["RANSAC"] = BenchmarkGroup()

function setupme(iterations)
    vs, ns, norms4Plot, shape_s = examplepc3()
    pcr = PointCloud(vs, ns, 32)
    p = RANSACParameters{Float64}()
    # plane
	p = RANSACParameters(p, ϵ_plane=0.3, α_plane=deg2rad(5))
    # cylidner
	p = RANSACParameters(p, ϵ_cylinder=0.3, α_cylinder=deg2rad(5))
    # sphere
	p = RANSACParameters(p, ϵ_sphere=0.3, α_sphere=deg2rad(5))
    # number of minimal subsets drawed in one iteration
	# probability that we found shapes
	# minimum shape size
	# maximum number of iteration
	p = RANSACParameters(p, minsubsetN=15, prob_det=0.9, τ=900)
	p = RANSACParameters(p, itermax=iterations, drawN=3)
    return (pcr, p, true)
end

Random.seed!(1234);
sharp = setupme(100_000)
SUITE["RANSAC"]["smallideal"] = @benchmarkable ransac($sharp..., reset_rand=true)
