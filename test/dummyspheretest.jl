# Dummy test of sphere fitting

using RANSAC: issphere, isplane
using RANSAC: FittedPlane, FittedSphere
using RANSAC: RANSACParameters

const EPSI = 0.1
const ALFI = deg2rad(10)
rp = RANSACParameters{Float64}()
const defrp = RANSACParameters(rp, ϵ_sphere = EPSI, α_sphere = ALFI)

@testset "true sphere 1" begin
    tv1 = [SVector(0,-1,0.0), SVector(0,0,-1.0), SVector(1,0,0.0), SVector(0,1,0.0)]
    tn1 = [SVector(0,-1,0.0), SVector(0,0,-1.0), SVector(1,0,0.0), SVector(0,1,0.0)]
    fs = issphere(tv1, tn1, defrp)
    fp = isplane(tv1, tn1, RANSACParameters(defrp, α_plane=π/2, collin_threshold=0.2))
    @test fs.issphere == true
    @test fp.isplane == false
end

@testset "true sphere 2" begin
    tv1 = [SVector(0,-0.99,0.0), SVector(0,0,-1.0), SVector(1.01,0,0.0), SVector(0,1,0.0)]
    tn1 = [SVector(0,-1,0.0), SVector(0,0,-1.0), SVector(1,0,0.0), SVector(0,1,0.0)]
    fs1 = issphere(tv1, tn1, defrp) # true
    fs2 = issphere(tv1, tn1, RANSACParameters(defrp, ϵ_sphere=0.01)) # false, cause ϵ
    fp = isplane(tv1, tn1, RANSACParameters(defrp, α_plane=π/2, collin_threshold=0.2))
    @test fs1.issphere == true
    @test fs2.issphere == false
    @test fp.isplane == false
end

@testset "false sphere 1" begin
    tv1 = [SVector(0,1,0.0), SVector(0,0,-1.0), SVector(1,0,0.0), SVector(0,1,0.0)]
    tn1 = [SVector(0,-1,0.0), SVector(0,0,-1.0), SVector(1,0,0.0), SVector(0,1,0.0)]
    fs1 = issphere(tv1, tn1, defrp)
    # should be false even with large thresholds
    fs2 = issphere(tv1, tn1, RANSACParameters(defrp, ϵ_sphere=10, α_sphere=π/2))
    fp = isplane(tv1, tn1, RANSACParameters(defrp, α_plane=π/2, collin_threshold=0.2))
    @test fs1.issphere == false
    @test fs2.issphere == false
    @test fp.isplane == false
end
