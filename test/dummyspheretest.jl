# Dummy test of sphere fitting

using RANSAC: FittedPlane, FittedSphere, fit
#using RANSAC: RANSACParameters

const EPSI = 0.1
const ALFI = deg2rad(10)
rp = ransacparameters()
#const defrp = RANSACParameters(rp, ϵ_sphere = EPSI, α_sphere = ALFI)
const defrp = ransacparameters(rp; sphere=(ϵ=EPSI, α=ALFI,))

@testset "true sphere 1" begin
    tv1 = [SVector(0,-1,0.0), SVector(0,0,-1.0), SVector(1,0,0.0), SVector(0,1,0.0)]
    tn1 = [SVector(0,-1,0.0), SVector(0,0,-1.0), SVector(1,0,0.0), SVector(0,1,0.0)]
    fs = fit(FittedSphere, tv1, tn1, defrp)
    #fp = fitplane(tv1, tn1, RANSACParameters(defrp, α_plane=π/2, collin_threshold=0.2))
    fp = fit(FittedPlane, tv1, tn1, ransacparameters(defrp, plane=(α=π/2,), common=(collin_threshold=0.2,)))
    @test fs isa FittedSphere
    @test fp === nothing
end

@testset "true sphere 2" begin
    tv1 = [SVector(0,-0.99,0.0), SVector(0,0,-1.0), SVector(1.01,0,0.0), SVector(0,1,0.0)]
    tn1 = [SVector(0,-1,0.0), SVector(0,0,-1.0), SVector(1,0,0.0), SVector(0,1,0.0)]
    fs1 = fit(FittedSphere, tv1, tn1, defrp) # true
    #fs2 = fitsphere(tv1, tn1, RANSACParameters(defrp, ϵ_sphere=0.01)) # false, cause ϵ
    fs2 = fit(FittedSphere, tv1, tn1, ransacparameters(defrp, sphere=(ϵ=0.01,))) # false, cause ϵ
    #fp = fitplane(tv1, tn1, RANSACParameters(defrp, α_plane=π/2, collin_threshold=0.2))
    fp = fit(FittedPlane, tv1, tn1, ransacparameters(defrp, plane=(α=π/2,), common=(collin_threshold=0.2,)))
    @test fs1 isa FittedSphere
    @test fs2 === nothing
    @test fp === nothing
end

@testset "false sphere 1" begin
    tv1 = [SVector(0,1,0.0), SVector(0,0,-1.0), SVector(1,0,0.0), SVector(0,1,0.0)]
    tn1 = [SVector(0,-1,0.0), SVector(0,0,-1.0), SVector(1,0,0.0), SVector(0,1,0.0)]
    fs1 = fit(FittedSphere, tv1, tn1, defrp)
    # should be false even with large thresholds
    #fs2 = fitsphere(tv1, tn1, RANSACParameters(defrp, ϵ_sphere=10, α_sphere=π/2))
    fs2 = fit(FittedSphere, tv1, tn1, ransacparameters(defrp, sphere=(ϵ=10,α=π/2,)))
    #fp = fitplane(tv1, tn1, RANSACParameters(defrp, α_plane=π/2, collin_threshold=0.2))
    fp = fit(FittedPlane, tv1, tn1, ransacparameters(defrp, plane=(α=π/2,), common=(collin_threshold=0.2,)))
    @test fs1 === nothing
    @test fs2 === nothing
    @test fp === nothing
end
