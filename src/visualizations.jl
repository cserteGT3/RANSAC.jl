"""
    showgeometry(scene, vs, ns; arrow = 0.5)

Show pointcloud with normals.
"""
function showgeometry(vs, ns; arrow = 0.5)
        plns = normalsforplot(vs, ns, arrow)
        scene = scatter(vs)
        linesegments!(scene, plns, color = :blue)
end

function showcandlength(ck)
    for c in ck
        println("candidate length: $(length(c.inpoints))")
    end
end

function showshapes(s, pointcloud, candidateA)
    colA = [:blue, :black, :darkred, :green, :brown, :yellow, :orange, :lightsalmon1, :goldenrod4, :olivedrab2, :indigo, :lightgreen, :darkorange1, :green2]
    @assert length(candidateA) <= length(colA) "Not enough color in colorarray. Fix it manually. :/"
    for i in 1:length(candidateA)
        ind = candidateA[i].inpoints
        scatter!(s, pointcloud.vertices[ind], color = colA[i])
    end
    s
end

function showshapes(pointcloud, candidateA)
    sc = Scene()
    showshapes(sc, pointcloud, candidateA)
end

function getrest(pc)
    return findall(pc.isenabled)
end


function showtype(l)
    for t in l
        println(t.candidate.shape)
    end
end

function showbytype(s, pointcloud, candidateA)
    for c in candidateA
        ind = c.inpoints
        if c.candidate.shape isa FittedCylinder
            colour = :red
        elseif c.candidate.shape isa FittedSphere
            colour = :green
        elseif c.candidate.shape isa FittedPlane
            colour = :orange
        end
        scatter!(s, pointcloud.vertices[ind], color = colour)
    end
    s
end

function showbytype(pointcloud, candidateA)
    sc = Scene()
    showbytype(sc, pointcloud, candidateA)
end
