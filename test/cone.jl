@testset "easy test" begin
    apex = SVector(15.0,-7, 0.1)
    axis = SVector(0,0,-1.)
    opang = deg2rad(20)
    h = 10.
    sizet = (10,10)

    cps, cns = samplecone(apex, axis, opang, h, sizet)

    rrr = [25, 19, 97]
    hicone = RANSAC.fit3pointcone(cps[rrr], cns[rrr])
    rrr2 = [25, 97, 19]
    hicone2 = RANSAC.fit3pointcone(cps[rrr2], cns[rrr2])

    @test hicone.apex-apex |> norm < 1e-13
    @test hicone2.apex-apex |> norm < 1e-13

    @test hicone.axis+[0,0,1] |> norm < 1e-13
    @test hicone2.axis+[0,0,1] |> norm < 1e-13

    @test isapprox(hicone.opang, opang)
    @test isapprox(hicone2.opang, opang)

    p = RANSACParameters{Float64}()
    hcone = RANSAC.fitcone(cps[rrr], cns[rrr], p)
    hcone2 = RANSAC.fitcone(cps[rrr2], cns[rrr2], p)

    @test hcone isa FittedCone
    @test hcone2 isa FittedCone
end
