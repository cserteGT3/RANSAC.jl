# test yaml config file import

## test isntequal()
@testset "test isntequal itself" begin
    nt1 = (a=1.5, b=3, c=[1, 2, 3])
    nt2 = (c=[1, 2, 3], b=3, a=1.5)
    @test isntequal(nt1, nt1)
    @test isntequal(nt1, nt2)

    nt3 = (c=[1, 2, 3], b=3.0, a=1.5)
    @test isntequal(nt1, nt3)
    @test isntequal(nt2, nt3)
    
    nt4 = (c=[1, 2, 3], d=3.1, a=1.5)
    nt5 = (c=[1, 2, 3], d=3, a=1.5)
    @test ! isntequal(nt1, nt4)
    @test ! isntequal(nt2, nt4)
    @test ! isntequal(nt3, nt4)
    @test ! isntequal(nt1, nt5)
    @test ! isntequal(nt2, nt5)
    @test ! isntequal(nt3, nt5)
    @test ! isntequal(nt4, nt5)

    # nested name tuples
    nt11 = (sphere = nt1, beta=2, gamma="str1")
    nt1_ = (beta=2, gamma="str1", sphere = nt1)
    nt13 = (beta=2, sphere = nt3, gamma="str1")
    nt12 = (sphere = nt2, gamma="str1", beta=2)
    eqnts = [nt11, nt1_, nt13, nt12]
    eqnts_size = size(eqnts, 1)
    for i in 1:eqnts_size
        for j in i:eqnts_size
            @test isntequal(eqnts[i], eqnts[j])
        end
    end

    # not equal nesteds
    nnt11 = (sphere = nt1, beta=2.1, gamma="str1")
    nnt1_ = (beta=2, gamma="str2", sphere = nt1)
    nnt13 = (beta=2, sphere = nt4, gamma="str1")
    nnt12 = (sphere = nt4, gamma="str1", beta=2.1)
    nnt15 = (sphere = nt5, name="str1", beta=2.1)

    noteqnts = [nnt11, nnt1_, nnt13, nnt12, nnt15]
    noteqnts_size = size(noteqnts, 1)
    for i in 1:eqnts_size
        for j in i:noteqnts_size
            @test ! isntequal(eqnts[i], noteqnts[j])
        end
    end

    # nest moar
    nt111 = (squirrel = nt11, beta=2, gamma="str1")
    nt133 = (beta=2, squirrel = nt13, gamma="str1")
    @test isntequal(nt111, nt133)
end

@testset "yaml config read - t1.yml" begin
    f1 = joinpath(pwd(), "yaml", "t1.yml")    
    
    conf1 = readconfig(f1)

    p1 = ransacparameters()
    p1 = ransacparameters(p1,; sphere=(ϵ=0.2, α=0.05, sphere_par=0.01,), plane=(ϵ=0.1, α=0.01,))
    p1 = ransacparameters(p1,; cylinder=(α=0.0872,), cone=(ϵ=1, α=3.14, minconeopang=1.,))
    p1 = ransacparameters(p1,; iteration=(drawN=9, minsubsetN=2, prob_det=0.999,τ=10000, itermax=100000, shape_types=[FittedPlane, FittedSphere],))
    p1 = ransacparameters(p1,; common=(parallelthrdeg=0.5, collin_threshold=0.3,))

    @test isntequal(conf1, p1)
end

@testset "yaml config read - t2.yml" begin
    f1 = joinpath(pwd(), "yaml", "t2.yml")    
    
    conf1 = readconfig(f1)

    p1 = ransacparameters()
    p1 = ransacparameters(p1, plane=(ϵ=0.35, α=1.0872,), sphere=(sphere_par=0.025,))
    p1 = ransacparameters(p1, iteration=(itermax=100,))
    p1 = ransacparameters(p1, common=(collin_threshold=0.22, parallelthrdeg=1.2,))

    @test isntequal(conf1, p1)
end