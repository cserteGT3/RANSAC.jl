@testset "IterationCandidates" begin
    ic = RANSAC.IterationCandidates()
    fp = FittedPlane(SVector(0.5,0.5,0.5), SVector{3}(0,0,1.))
    sc = ConfidenceInterval(0,1)
    inds = [1,2,3,4,5]
    
    @test length(ic) == 0
    RANSAC.recordscore!(ic, fp, sc, inds)
    @test length(ic) == 1
    @test size(ic.shapes,1) == 1
    @test size(ic.scores,1) == 1
    @test size(ic.inpoints,1) == 1

    deleteat!(ic, 1)
    @test length(ic) == 0
    @test size(ic.shapes,1) == 0
    @test size(ic.scores,1) == 0
    @test size(ic.inpoints,1) == 0
end