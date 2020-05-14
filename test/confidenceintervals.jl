@testset "confident tests" begin
    a = 1.0
    b = 3
    ci1 = ConfidenceInterval(a, b)
    @test ci1.E == 2.0
    @test ci1.min === 1.0
    @test ci1.max === 3.0

    @test_throws ErrorException ConfidenceInterval(b, a)
end

@testset "not so confident tests" begin
    a = 9.7
    b = 153.9

    nc1 = notsoconfident(b, a)
    nc2 = notsoconfident(a, b)
    
    @test nc1.min == 9.7
    @test nc1.max == 153.9
    @test nc1.E == 81.8
        
    @test nc1.min == nc2.min
    @test nc1.max == nc2.max
    @test nc1.E == nc2.E
end