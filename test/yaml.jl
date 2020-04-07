# test yaml config file import

function compareallfields(p1, p2)
    @test typeof(p1) == typeof(p2)
    fnames = fieldnames(typeof(p1))
    for fn in fnames
        @test getproperty(p1, fn) == getproperty(p2, fn)
    end
end

@testset "compareallfields itself" begin
    p1 = RANSACParameters()
    compareallfields(p1, p1)
end

@testset "yaml config read - all fields specified" begin
    f1 = joinpath(pwd(), "yaml", "t1.yml")    
    
    conf1 = readconfig(f1)
    confmk1 = RANSACParameters( ϵ_plane=0.1,
                                α_plane=0.01,
                                ϵ_sphere=0.2,
                                α_sphere=0.05,
                                ϵ_cylinder=0.3,
                                α_cylinder=0.0872,
                                ϵ_cone=1.,
                                α_cone=3.14,
                                minconeopang=1.,
                                ϵ_torus=0.3,
                                α_torus=0.0872,
                                drawN=5,
                                minsubsetN=2,
                                prob_det=0.99,
                                τ=10000,
                                itermax=100_000,
                                parallelthrdeg=0.5,
                                collin_threshold=0.3,
                                β=1.,
                                sphere_par=0.01,
                                shape_types=[:plane])
    compareallfields(conf1, confmk1)

    # Float32
    confmk1_f32 = RANSACParameters{Float32}( ϵ_plane=0.1,
                                α_plane=0.01,
                                ϵ_sphere=0.2,
                                α_sphere=0.05,
                                ϵ_cylinder=0.3,
                                α_cylinder=0.0872,
                                ϵ_cone=1.,
                                α_cone=3.14,
                                minconeopang=1.,
                                ϵ_torus=0.3,
                                α_torus=0.0872,
                                drawN=5,
                                minsubsetN=2,
                                prob_det=0.99,
                                τ=10000,
                                itermax=100_000,
                                parallelthrdeg=0.5,
                                collin_threshold=0.3,
                                β=1.,
                                sphere_par=0.01,
                                shape_types=[:plane])
    conf1_f32 = readconfig(f1, RANSACParameters{Float32})
    compareallfields(conf1_f32, confmk1_f32)


end

@testset "yaml config read - some fields specified" begin
    f2 = joinpath(pwd(), "yaml", "t2.yml")    
    
    conf2 = readconfig(f2)
    confmk2 = RANSACParameters( ϵ_plane=0.35,
                                α_plane=1.0872,
                                itermax=100,
                                parallelthrdeg=1.2,
                                collin_threshold=0.22,
                                β=1.001,
                                sphere_par=0.025)
    compareallfields(conf2, confmk2)

    # Float32
    confmk2_f32 = RANSACParameters{Float32}( ϵ_plane=0.35,
                                α_plane=1.0872,
                                itermax=100,
                                parallelthrdeg=1.2,
                                collin_threshold=0.22,
                                β=1.001,
                                sphere_par=0.025)
    conf2_f32 = readconfig(f2, RANSACParameters{Float32})
    compareallfields(conf2_f32, confmk2_f32)
end