# fitting

"""
    struct FittedCylinder{A<:AbstractArray, R<:Real} <: FittedShape

Cylinder primitive, defined by its axis direction,
a point that lies on its axis,
and its radius.
Also stored, if the normals point outwards of the shape.
"""
struct FittedCylinder{A<:AbstractArray, R<:Real} <: FittedShape
    axis::A
    center::A
    radius::R
    outwards::Bool
end

Base.show(io::IO, x::FittedCylinder) =
    print(io, """cylinder, R: $(x.radius)""")

Base.show(io::IO, ::MIME"text/plain", x::FittedCylinder{A, R}) where {A, R} =
    print(io, """FittedCylinder{$A, $R}\ncylinder, center: $(x.center), axis: $(x.axis), R: $(x.radius), $(x.outwards ? "outwards" : "inwards" )""")

strt(x::FittedCylinder) = "cylinder"

function setcylinderOuterity(fc, b)
    FittedCylinder(fc.axis, fc.center, fc.radius, b)
end

function fit2pointcylinder(p, n, params)
    @unpack α_cylinder, ϵ_cylinder, parallelthrdeg = params
    # the normals are too parallel so cannot fit cylinder to it
    if abs(dot(n[1], n[2])) > cosd(parallelthrdeg)
        return nothing
    end

    an = normalize(n[1] × n[2])

    function lineplaneintersect(n, u, w)
        # http://geomalgorithms.com/a05-_intersect-1.html
        # n=plane normal
        # u= direction of the line
        # w=

        return dot(-n, w)/dot(n,u)
    end

    function project2plane(n, w)
        #n: normal
        #w: helyvektor
        return w+n*lineplaneintersect(n,n,w)
    end

    function projectto2d(xaxis, yaxis, zaxis, p1)
        # xax: pl a ponban muataó vektor
        xx = xaxis[1];
        xy = xaxis[2];
        xz = xaxis[3];

        yx = yaxis[1];
        yy = yaxis[2];
        yz = yaxis[3];

        zx = zaxis[1];
        zy = zaxis[2];
        zz = zaxis[3];

        px = p1[1];
        py = p1[2];
        pz = p1[3];


        r1 = -((-(pz*yy*zx) + py*yz*zx + pz*yx*zy - px*yz*zy - py*yx*zz + px*yy*zz)/(xz*yy*zx - xy*yz*zx - xz*yx*zy + xx*yz*zy + xy*yx*zz - xx*yy*zz))
        r2 = -((pz*xy*zx - py*xz*zx - pz*xx*zy + px*xz*zy + py*xx*zz - px*xy*zz)/(xz*yy*zx - xy*yz*zx - xz*yx*zy + xx*yz*zy + xy*yx*zz - xx*yy*zz))
        r3 = -((pz*xy*yx - py*xz*yx - pz*xx*yy + px*xz*yy + py*xx*yz - px*xy*yz)/(-(xz*yy*zx) + xy*yz*zx + xz*yx*zy - xx*yz*zy - xy*yx*zz + xx*yy*zz))

        return SVector(r1, r2)
    end

    function lineintersectionpoint(l1, l2)
        # l1 = [[p1x, p1y],[p2x,p2y]]
        a = l1[1]
        b = l1[2]
        c = l2[1]
        d = l2[2]

        amb_ = a-b
        cmd_ = c-d

        d1 = det(vcat(a',b'))
        d2 = det(vcat(c',d'))
        d3 = det(vcat(amb_',cmd_'))
        return (d1*cmd_-d2*amb_)/d3
    end

    xax = normalize(project2plane(an, p[1]))
    yax = normalize(an × xax)

    p11proj = projectto2d(xax, yax, an, project2plane(an, p[1]))
    p12proj = projectto2d(xax, yax, an, project2plane(an, p[1]+n[1]))

    p21proj = projectto2d(xax, yax, an, project2plane(an, p[2]))
    p22proj = projectto2d(xax, yax, an, project2plane(an, p[2]+n[2]))

    interc = lineintersectionpoint([p11proj,p12proj], [p21proj,p22proj])

    c = interc[1]*xax+interc[2]*yax

    nnormies = [norm(pt - c - an*dot( an, pt-c )) for pt in p[1:2]]
    R = (nnormies[1] + nnormies[2])/2
    # R = (norm(interc-p11proj) + norm(interc-p21proj))/2

    #p[1] - an*dot( an, p[1]-c )

    outw = dot(p12proj-p11proj, p11proj-interc) > 0

    return FittedCylinder(an, c, R, outw)
end

"""
    fitcylinder(p, n, params)

Fit a cylinder to 2 points. Additional points and their normals are used to validate the fit.
Normals are expected to be normalized.
Return `nothing` if points do not fit to a cylinder.
"""
function fitcylinder(p, n, params)
    @unpack ϵ_cylinder, α_cylinder, parallelthrdeg = params
    pl = length(p)
    @assert pl == length(n) "Size must be the same."
    @assert pl > 2 "Size must be at least 3."

    # "forcefit" a cylinder
    fc = fit2pointcylinder(p, n, params)

    fc === nothing && return nothing


    thr = cos(α_cylinder)
    vert_ok = falses(pl)
    norm_ok = falses(pl)
    invnorm_ok = falses(pl)
    for i in 1:pl
        # current normal
        curr_norm = p[i] - fc.axis*dot( fc.axis, p[i]-fc.center ) - fc.center
        # vertice check
        vert_ok[i] = abs(norm(curr_norm)-fc.radius) < ϵ_cylinder
        # normal check
        dotp = dot( normalize(curr_norm), n[i] )
        norm_ok[i] = dotp > thr
        invnorm_ok[i] = dotp < -thr
    end
    vert_ok == trues(pl) || return nothing
    norm_ok == trues(pl) && return setcylinderOuterity(fc, true)
    invnorm_ok == trues(pl) && return setcylinderOuterity(fc, false)
    return nothing
end

# bitmapping

function scorecandidate(pc, candidate::ShapeCandidate{T}, subsetID, params) where {T<:FittedCylinder}
    ps = @view pc.vertices[pc.subsets[subsetID]]
    ns = @view pc.normals[pc.subsets[subsetID]]
    ens = @view pc.isenabled[pc.subsets[subsetID]]

    cp, pp = compatiblesCylinder(candidate.shape, ps, ns, params)
    inder = cp.&ens
    inpoints = (pc.subsets[subsetID])[inder]
    #inpoints = ((pc.subsets[1])[ens])[cp]
    score = estimatescore(length(pc.subsets[subsetID]), pc.size, length(inpoints))
    pc.levelscore[candidate.octree_lev] += E(score)
    return ScoredShape(candidate, score, inpoints)
end

"""
    compatiblesCylinder(cylinder, points, normals, eps, alpharad)

Create a bool-indexer array for those points that are compatible to the cylinder.
Give back the projected points too for parameter space magic.

Compatibility is measured with an `eps` distance to the cylinder and an `alpharad` angle to it's normal.
"""
function compatiblesCylinder(cylinder, points, normals, params)
    @unpack ϵ_cylinder, α_cylinder = params
    @assert length(points) == length(normals) "Size must be the same."

    c = cylinder.center
    R = cylinder.radius
    a = cylinder.axis

    comp = falses(length(points))
    pars = fill(SVector(0.0,0.0), length(points))

    for i in 1:length(points)
        curr_norm = points[i] - a*dot( a, points[i]-c ) - c
        # if the radius is correct
        if abs(norm(curr_norm)-R) < ϵ_cylinder
            if cylinder.outwards
                comp[i] = isparallel(normalize(curr_norm), normals[i], α_cylinder)
            else
                comp[i] = isparallel(-normalize(curr_norm), normals[i], α_cylinder)
            end
        end

        # playing with parameters:
        #TODO: implement it
    end
    return comp, pars
end

"""
    refit(s, pc, ϵ, α)

Refit cylinder. Only s.inpoints is updated.
"""
function refitcylinder(s, pc, params)
    # TODO: use octree for that
    cp, _ = compatiblesCylinder(s.candidate.shape, pc.vertices[pc.isenabled], pc.normals[pc.isenabled], params)
    s.inpoints = ((1:pc.size)[pc.isenabled])[cp]
    s
end
