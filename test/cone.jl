@testset "easy test" begin
    apex = SVector(0.,0,0)
    axis = SVector(0,0,-1.)
    opang = deg2rad(20)
    h = 10.
    sizet = (10,10)

    p = RANSACParameters{Float64}()

    cps, cns = samplecone(apex, axis, opang, h, sizet)

    rrr = [25, 19, 97]
    hicone = RANSAC.fit3pointcone(cps[rrr], cns[rrr], p)
    rrr2 = [25, 97, 19]
    hicone2 = RANSAC.fit3pointcone(cps[rrr2], cns[rrr2], p)

    @test hicone.apex |> norm < 1e-14
    @test hicone2.apex |> norm < 1e-14

    @test hicone.axis+[0,0,1] |> norm < 1e-14
    @test hicone2.axis+[0,0,1] |> norm < 1e-14

    @test isapprox(hicone.opang, opang)
    @test isapprox(hicone2.opang, opang)
end
