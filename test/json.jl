# testing json export by
# exporting, then importing
@testset "all in one" begin
    # toDict

    ## FottedShape
    
    # plane
    p1 = SVector(15.6, 0, -13.7)
    n1 = normalize(SVector(34,45,7))
    s_plane = FittedPlane(p1, n1)
    d_plane = RANSAC.toDict(s_plane)
    
    @test d_plane == Dict("type"=>"plane", "point"=>p1, "normal"=>n1)

    # sphere
    s_sphere = FittedSphere(p1, 13.23444, true)
    d_sphere = RANSAC.toDict(s_sphere)

    @test d_sphere == Dict("type"=>"sphere", "radius"=>13.23444, "center"=>p1, "outwards"=>true)

    # cylinder
    a1 = SVector(-17.1, 8, 2.42)
    s_cylinder = FittedCylinder(a1, p1, 0.13, false)
    d_cylinder = RANSAC.toDict(s_cylinder)

    @test d_cylinder == Dict("type"=>"cylinder", "axis"=>a1, "center"=>p1, "radius"=>0.13, "outwards"=>false)

    # cone
    ap1 = SVector(0.0,0,0)
    ax1 = normalize(SVector(-1.5, 7, 2))
    s_cone = FittedCone(ap1, ax1, 0.785, true)
    d_cone = RANSAC.toDict(s_cone)

    @test d_cone == Dict("type"=>"cone", "apex"=>ap1, "axis"=>ax1, "opang"=>0.785, "outwards"=>true)

    ## ShapeCandidate
    sc1 = ShapeCandidate(s_plane, 5)
    @test RANSAC.toDict(sc1) == d_plane

    ## ScoredShape
    ss1 = ScoredShape(ShapeCandidate(s_cone, 0), ConfidenceInterval(0,1.0), [1,2,3])
    @test RANSAC.toDict(ss1) == d_cone

    ## Array of primitives
    sa = [s_plane, s_sphere, s_cylinder, s_cone]
    sa_dict = [d_plane, d_sphere, d_cylinder, d_cone]
    
    d_sa = RANSAC.toDict(sa)

    @test d_sa == Dict("primitives"=>sa_dict)

    ## Array of ShapeCandidate
    sc_a = [sc1]
    sc_a_dict = [d_plane]

    @test RANSAC.toDict(sc_a) == Dict("primitives"=>sc_a_dict)

    ## Array of ScoredShape
    ss2 = ScoredShape(ShapeCandidate(s_cylinder, 0), ConfidenceInterval(0.123,0.5), [1,2,3, 4])
    ss_a = [ss1, ss2]
    ss_a_dict = [d_cone, d_cylinder]

    @test RANSAC.toDict(ss_a) == Dict("primitives"=>ss_a_dict)


    # exportJSON

    ## with indentation
    io1 = IOBuffer()
    exportJSON(io1, s_cone, 2)
    string_io1 = String(take!(io1))
    #@test string_io1 == "{\n  \"opang\": 0.785,\n  \"outwards\": true,\n  \"axis\": [\n    -0.20180183819889375,\n    0.9417419115948374,\n    0.269069117598525\n  ],\n  \"apex\": [\n    0.0,\n    0.0,\n    0.0\n  ],\n  \"type\": \"cone\"\n}\n"
    pd1 = JSON.parse(string_io1)
    @test pd1 == d_cone

    io2 = IOBuffer()
    exportJSON(io2, [s_plane, s_cone], 1)
    string_io2 = String(take!(io2))
    #@test string_io2 == "{\n \"primitives\": [\n  {\n   \"point\": [\n    15.6,\n    0.0,\n    -13.7\n   ],\n   \"normal\": [\n    0.5982430416161189,\n    0.7917922609625102,\n    0.1231676850386127\n   ],\n   \"type\": \"plane\"\n  },\n  {\n   \"opang\": 0.785,\n   \"outwards\": true,\n   \"axis\": [\n    -0.20180183819889375,\n    0.9417419115948374,\n    0.269069117598525\n   ],\n   \"apex\": [\n    0.0,\n    0.0,\n    0.0\n   ],\n   \"type\": \"cone\"\n  }\n ]\n}\n"
    pd2 = JSON.parse(string_io2)
    @test pd2 == Dict("primitives"=>[d_plane, d_cone])

    ## Without indentation
    io3 = IOBuffer()
    exportJSON(io3, s_cylinder)
    string_io3 = String(take!(io3))
    #@test string_io3 == "{\"outwards\":false,\"axis\":[-17.1,8.0,2.42],\"radius\":0.13,\"center\":[15.6,0.0,-13.7],\"type\":\"cylinder\"}"
    pd3 = JSON.parse(string_io3)
    @test pd3 == d_cylinder

    io4 = IOBuffer()
    exportJSON(io4, [s_plane, s_sphere, s_cone, s_plane])
    string_io4 = String(take!(io4))
    #@test string_io4 == "{\"primitives\":[{\"point\":[15.6,0.0,-13.7],\"normal\":[0.5982430416161189,0.7917922609625102,0.1231676850386127],\"type\":\"plane\"},{\"outwards\":true,\"radius\":13.23444,\"center\":[15.6,0.0,-13.7],\"type\":\"sphere\"},{\"opang\":0.785,\"outwards\":true,\"axis\":[-0.20180183819889375,0.9417419115948374,0.269069117598525],\"apex\":[0.0,0.0,0.0],\"type\":\"cone\"},{\"point\":[15.6,0.0,-13.7],\"normal\":[0.5982430416161189,0.7917922609625102,0.1231676850386127],\"type\":\"plane\"}]}"
    pd4 = JSON.parse(string_io4)
    @test pd4 == Dict("primitives"=>[d_plane, d_sphere, d_cone, d_plane])

    ## Write to file w&wo indentation
    f1, fio1 = mktemp(pwd())
    close(fio1)
    open(f1, "w") do io
        exportJSON(io, [s_plane, s_sphere, s_cone, s_plane], 4)
    end
    str_f1 = read(f1, String)
    pd_f1 = JSON.parse(str_f1)
    @test pd_f1 == Dict("primitives"=>[d_plane, d_sphere, d_cone, d_plane])
    Base.Filesystem.rm(f1)

    f2, fio2 = mktemp(pwd())
    close(fio2)
    open(f2, "w") do io
        exportJSON(io, [s_plane, s_sphere, s_cone, s_plane, s_cylinder])
    end
    str_f2 = read(f2, String)
    pd_f2 = JSON.parse(str_f2)
    @test pd_f2 == Dict("primitives"=>[d_plane, d_sphere, d_cone, d_plane, d_cylinder])
    Base.Filesystem.rm(f2)
end