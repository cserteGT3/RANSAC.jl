@testset "segmentdistance" begin
    segm = [[2.0,2.0], [4.0,2.0]]

    p1 = [0.0,2.0];
    @test isapprox(RANSAC.segmentdistance(p1, segm), 2)
    p1 = [2.0,0.0];
    @test isapprox(RANSAC.segmentdistance(p1, segm), 2)
    p1 = [2.0,2.0];
    @test isapprox(RANSAC.segmentdistance(p1, segm), 0)
    p1 = [3.0,2.0];
    @test isapprox(RANSAC.segmentdistance(p1, segm), 0)
    p1 = [4.0,2.0];
    @test isapprox(RANSAC.segmentdistance(p1, segm), 0)
    p1 = [5.0,2.0];
    @test isapprox(RANSAC.segmentdistance(p1, segm), 1)

    p1 = [3,1.0];
    @test isapprox(RANSAC.segmentdistance(p1, segm), 1)
    p1 = [3,-1.0];
    @test isapprox(RANSAC.segmentdistance(p1, segm), 3)

    p1 = [3,3];
    @test isapprox(RANSAC.segmentdistance(p1, segm), 1)
end
