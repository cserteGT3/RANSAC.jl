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

function plotshape(shape::FittedShape; kwargs...)
    plotshape!(Scene(), shape; kwargs...)
end

function plotshape!(sc, shape::FittedPlane; scale=(1.,1.), color=(:blue, 0.1))
    # see project2plane
    o_z = normalize(shape.normal)
    o_x = normalize(arbitrary_orthogonal(o_z))
    o_y = normalize(cross(o_z, o_x))

    p1 = shape.point
    p2 = p1 + scale[1]*o_x
    p3 = p1 + scale[1]*o_x + scale[2]*o_y
    p4 = p1 + scale[2]*o_y

    mesh!(sc, [p1,p2,p3], color=color, transparency=true)
    mesh!(sc, [p1,p3,p4], color=color, transparency=true)
    sc
end

function plotshape!(sc, shape::FittedSphere; scale=(1.,), color=(:blue, 0.1))
    mesh!(sc, Sphere(Point(shape.center), scale[1]*shape.radius), color=color, transparency=true)
end

function plotshape!(sc, shape::FittedCylinder; scale=(1.,), color=(:blue, 0.1))
    o = Point(shape.center)
    extr = o+scale[1]*Point(normalize(shape.axis))
    mesh!(sc, Cylinder(o, extr, shape.radius), color=color, transparency=true)
end

function shiftplane!(sc, p::FittedPlane, dist; kwargs...)
    newo = p.point+dist*normalize(p.normal)
    newp = FittedPlane(true, newo, p.normal)
    plotshape!(sc, newp; kwargs...)
end
