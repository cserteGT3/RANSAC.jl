# Test iswithinrectangle function

using RANSAC: iswithinrectangle

using RegionTrees: HyperRectangle, vertices

rectangle = HyperRectangle(SVector(0.0,0,0), SVector(1.0,1,1))

@testset "corner points" begin
    # false
    v000 = SVector(0.0,0,0)
    v100 = SVector(1.0,0,0)
    v010 = SVector(0.0,1,0)
    v001 = SVector(0.0,0,1)
    v110 = SVector(1.0,1,0)
    v101 = SVector(1.0,0,1)
    v011 = SVector(0.0,1,1)

    # true
    v111 = SVector(1.0,1,1)

    @test iswithinrectangle(rectangle, v000) == false
    @test iswithinrectangle(rectangle, v100) == false
    @test iswithinrectangle(rectangle, v010) == false
    @test iswithinrectangle(rectangle, v001) == false
    @test iswithinrectangle(rectangle, v110) == false
    @test iswithinrectangle(rectangle, v101) == false
    @test iswithinrectangle(rectangle, v011) == false

    @test iswithinrectangle(rectangle, v111) == true
end

@testset "edge midpoints" begin
    # false
    v500 = SVector(0.5,0,0)
    v501 = SVector(0.5,0,1)
    v050 = SVector(0,0.5,0)
    v051 = SVector(0,0.5,1)
    v150 = SVector(1,0.5,0)
    v510 = SVector(0.5,1,0)
    v105 = SVector(1,0,0.5)
    v005 = SVector(0,0,0.5)
    v015 = SVector(0,1,0.5)

    # true
    v151 = SVector(1,0.5,1)
    v511 = SVector(0.5,1,1)
    v115 = SVector(1,1,0.5)

    @test iswithinrectangle(rectangle, v500) == false
    @test iswithinrectangle(rectangle, v501) == false
    @test iswithinrectangle(rectangle, v050) == false
    @test iswithinrectangle(rectangle, v051) == false
    @test iswithinrectangle(rectangle, v150) == false
    @test iswithinrectangle(rectangle, v510) == false
    @test iswithinrectangle(rectangle, v105) == false
    @test iswithinrectangle(rectangle, v005) == false
    @test iswithinrectangle(rectangle, v015) == false

    @test iswithinrectangle(rectangle, v151) == true
    @test iswithinrectangle(rectangle, v511) == true
    @test iswithinrectangle(rectangle, v115) == true
end

@testset "face midpoints" begin
    # false
    v550 = SVector(0.5,0.5,0)
    v505 = SVector(0.5,0,0.5)
    v055 = SVector(0,0.5,0.5)

    # true
    v551 = SVector(0.5,0.5,1)
    v515 = SVector(0.5,1,0.5)
    v155 = SVector(1,0.5,0.5)

    @test iswithinrectangle(rectangle, v550) == false
    @test iswithinrectangle(rectangle, v505) == false
    @test iswithinrectangle(rectangle, v055) == false

    @test iswithinrectangle(rectangle, v551) == true
    @test iswithinrectangle(rectangle, v515) == true
    @test iswithinrectangle(rectangle, v155) == true
end

@testset "inside points" begin
    # like facepoints but inside
    v5501 = SVector(0.5,0.5,0.1)
    v5015 = SVector(0.5,0.1,0.5)
    v0155 = SVector(0.1,0.5,0.5)
    v5509 = SVector(0.5,0.5,0.9)
    v5095 = SVector(0.5,0.9,0.5)
    v0955 = SVector(0.9,0.5,0.5)

    vin = [v5501, v5015, v0155, v5509, v5095, v0955]
    for vi in vin
        @test iswithinrectangle(rectangle, vi)
    end
end

@testset "outside points" begin
    # like facepoints but outside
    v5501m = SVector(0.5,0.5,-0.1)
    v5015m = SVector(0.5,-0.1,0.5)
    v0155m = SVector(-0.1,0.5,0.5)
    v5509p = SVector(0.5,0.5,1.1)
    v5095p = SVector(0.5,1.1,0.5)
    v0955p = SVector(1.1,0.5,0.5)

    vout = [v5501m, v5015m, v0155m, v5509p, v5095p, v0955p]
    for vo in vout
        @test iswithinrectangle(rectangle, vo) == false
    end
end
