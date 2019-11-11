@testset "largestconncomp on dense bitmap" begin
    pic1 = falses(150,150);
    pic1[55:75, 55:75] = trues(21,21);
    pic1[100:125, 100:125] = trues(26,26);
    pic1[130:140,130:140] = trues(11,11);
    indmap = [Int[] for i in 1:150, j in 1:150];
    for i in 55:75, j in 55:75
        append!(indmap[i,j], [1,2,3])
    end
    for i in 100:125, j in 100:125
        append!(indmap[i,j], [-99,-98])
    end
    for i in 130:140, j in 130:140
        append!(indmap[i,j], [0])
    end

    t1 = RANSAC.largestconncomp(pic1, indmap, 1:2)
    @test t1 == repeat([-99, -98], 26*26)
    @test t1 == RANSAC.largestconncomp(pic1, indmap)
    @test t1 == RANSAC.largestconncomp(pic1, indmap, trues(3,3))
end

@testset "8 connectivity" begin
    # largest patch is only 8 connectivity connected
    pic2 = falses(150,150);
    pic2[55:73, 55:73] = trues(19,19);
    pic2[100:2:126, 100:2:126] = trues(14,14);
    pic2[101:2:127, 101:2:127] = trues(14,14);
    pic2[130:140,130:140] = trues(11,11);
    indmap = [Int[] for i in 1:150, j in 1:150];
    for i in 55:73, j in 55:73
        append!(indmap[i,j], [1,2,3])
    end
    for i in 100:2:126, j in 100:2:126
        append!(indmap[i,j], [-99,-98])
    end
    for i in 101:2:127, j in 101:2:127
        append!(indmap[i,j], [-99,-98])
    end
    for i in 130:140, j in 130:140
        append!(indmap[i,j], [0])
    end

    t2d = RANSAC.largestconncomp(pic2, indmap)
    t2def = RANSAC.largestconncomp(pic2, indmap, :default)
    t24 = RANSAC.largestconncomp(pic2, indmap, 1:2)
    @test t24 == repeat([1,2,3], 19*19)
    @test t2d == t24
    @test t2def == t24

    t28 = RANSAC.largestconncomp(pic2, indmap, trues(3,3))
    t28eight = RANSAC.largestconncomp(pic2, indmap, :eight)
    @test t28 == repeat([-99,-98], 14*14*2)
    @test t28 == t28eight
end
