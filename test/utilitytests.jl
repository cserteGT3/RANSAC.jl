# Test functions from the Utilities module

using RANSAC: findAABB, smallestdistance

@testset "rand 3D points" begin
    rs3 = [SVector{3}(rand(3)) for _ in 1:30]
    iV3 = SVector(-1.0,-1,-1)
    aV3 = SVector(2.0,2.0,2.0)
    push!(rs3, iV3)
    push!(rs3, aV3)

    minV3, maxV3 = findAABB(rs3)
    @test minV3 == iV3
    @test maxV3 == aV3
end

@testset "rand 2D points" begin
    rs2 = [SVector{2}(rand(2)) for _ in 1:30]
    iV2 = SVector(-1.0,-1)
    aV2 = SVector(2.0,2.0)
    push!(rs2, iV2)
    push!(rs2, aV2)

    minV2, maxV2 = findAABB(rs2)
    @test minV2 == iV2
    @test maxV2 == aV2
end

@testset "smallest distance 2D" begin
    tv = [SVector(0.0, 0), SVector(1.0, 1), SVector(2.2, 2)]
    sd = smallestdistance(tv)
    @test isapprox(sd, sqrt(2))
end

@testset "smallest distance 3D" begin
    tv = [SVector(0.0,0,0), SVector(1.0,1,1), SVector(2.2,2,2)]
    sd = smallestdistance(tv)
    @test isapprox(sd, sqrt(3))
end

@testset "default common parameters" begin
    cp = (common=(collin_threshold = 0.2, parallelthrdeg = 1.),)
    @test isntequal(cp, RANSAC.defaultcommonparameters())
end

@testset "default iteration parameters" begin
    ip = (drawN=3, minsubsetN=15, prob_det=0.9, shape_types=[FittedSphere, FittedPlane, FittedCone], τ=900, itermax=1000, extract_s=:nofminset, terminate_s=:nofminset)
    @test isntequal((iteration=ip,), RANSAC.defaultiterationparameters([FittedSphere, FittedPlane, FittedCone]))

    ip = (drawN=3, minsubsetN=15, prob_det=0.9, shape_types=[FittedCone], τ=900, itermax=1000, extract_s=:nofminset, terminate_s=:nofminset)
    @test isntequal((iteration=ip,), RANSAC.defaultiterationparameters([FittedCone]))
end

@testset "default shape parameters" begin
    # sphere
    fsp = (sphere=(ϵ=0.3, α=deg2rad(5), sphere_par=0.02,),)
    @test isntequal(fsp, RANSAC.defaultshapeparameters(FittedSphere))

    # plane
    @test isntequal((plane=(ϵ=0.3, α=deg2rad(5)),), RANSAC.defaultshapeparameters(FittedPlane))
    # cylinder
    @test isntequal((cylinder=(ϵ=0.3, α=deg2rad(5)),), RANSAC.defaultshapeparameters(FittedCylinder))
    # cone
    @test isntequal((cone=(ϵ=0.3, α=deg2rad(5), minconeopang=deg2rad(2)),), RANSAC.defaultshapeparameters(FittedCone))
end

@testset "defaultparameters()" begin
    dfp = RANSAC.defaultparameters([FittedSphere])

    ip = (iteration=(drawN=3, minsubsetN=15, prob_det=0.9, shape_types=[FittedSphere], τ=900, itermax=1000, extract_s=:nofminset, terminate_s=:nofminset),)
    fsp = (sphere=(ϵ=0.3, α=deg2rad(5), sphere_par=0.02,),)
    cp = (common=(collin_threshold = 0.2, parallelthrdeg = 1.),)
    dip = merge(merge(ip, fsp), cp)
    @test isntequal(dfp, dip)

    dfp = RANSAC.defaultparameters([FittedSphere, FittedPlane])

    ip = (iteration=(drawN=3, minsubsetN=15, prob_det=0.9, shape_types=[FittedSphere, FittedPlane], τ=900, itermax=1000, extract_s=:nofminset, terminate_s=:nofminset),)
    dip = merge(dip, ip)
    dip = merge(dip, (plane=(ϵ=0.3, α=deg2rad(5)),))
    @test isntequal(dfp, dip)
end

@testset "ransacparameters() - named tuple" begin
    # named tuple version method
    rp = ransacparameters()
    @test isntequal(rp, RANSAC.DEFAULT_PARAMETERS)

    rp = ransacparameters(; sphere=(ϵ=0.9, α=deg2rad(1),), plane=(ϵ=1.0,))
    dp = RANSAC.defaultparameters([FittedPlane, FittedCone, FittedCylinder, FittedSphere])
    dp = merge(dp, (sphere=(ϵ=0.9, α=deg2rad(1),sphere_par = 0.02,),))
    dp = merge(dp, (plane=(ϵ=1.0,α=deg2rad(5),),))
    @test isntequal(rp, dp)
end

@testset "ransacparameters() - array" begin
    # method that takes an array of FittedShapes
    rp = ransacparameters([FittedCone, FittedCylinder])
    dpars = RANSAC.defaultparameters([FittedCone, FittedCylinder])

    @test isntequal(rp, dpars)

    rpn = ransacparameters([FittedCone, FittedCylinder], cone=(ϵ=0.9, α=deg2rad(12),), cylinder=(ϵ=0.001,))
    cop = RANSAC.defaultshapeparameters(FittedCone)
    copb = merge(cop.cone, (ϵ=0.9, α=deg2rad(12),))
    cop = (cone=copb,)

    cyp = RANSAC.defaultshapeparameters(FittedCylinder)
    cypb = merge(cyp.cylinder, (ϵ=0.001,))
    cyp = (cylinder=cypb,)
    dpars = merge(dpars, cyp)
    dpars = merge(dpars, cop)
    @test isntequal(rpn, dpars)
end

@testset "cross prod" begin
    A = rand(3,3)
    val = deg2rad(15)
    vec = normalize(rand(3))

    orival = A + sin(val) .* RANSAC.crossprodtensor(vec)
    dA = deepcopy(A)
    RANSAC.pluscrossprod!(dA, sin(val), vec)
    @test orival == dA

    orival2 = A + RANSAC.crossprodtensor(vec)
    dA2 = deepcopy(A)
    RANSAC.pluscrossprod!(dA2, 1, vec)
    @test orival2 == dA2

    newval3 = RANSAC.pluscrossprod!(deepcopy(A), 0, vec)
    @test A == newval3
end

@testset "push2candidatesandlevels!" begin
    fp = FittedPlane(SVector(0.5,0.5,0.5), SVector{3}(0,0,1.))
    candidates = FittedShape[]
    levels = Int[]

    RANSAC.push2candidatesandlevels!(candidates, fp, levels, 3)
    @test size(candidates) == (1,)
    @test size(levels) == (1,)

    RANSAC.push2candidatesandlevels!(candidates, [fp, fp], levels, 0)
    @test size(candidates) == (3,)
    @test size(levels) == (3,)
    @test levels[1] == 3
    @test levels[2] == 0
    @test levels[3] == 0
end