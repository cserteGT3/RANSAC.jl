"""
    project2plane(plane, points)

Project `points` on to the `plane`.

# Arguments:
- 'plane::FittedPlane': a plane.
- 'points::AbstractArray': an array of points (array-of-arrays).
"""
function project2plane(plane, points)
    # create a coordinate frame
    # z is the plane's normal
    o_z = normalize(plane.normal)
    # x is a random orthogonal vector (in the plane)
    o_x = normalize(arbitrary_orthogonal(o_z))
    # y is created so it's a right hand coord. frame
    o_y = normalize(cross(o_z, o_x))
    answer = similar(points)
    # get the coordinates of the points in the prev. created coord. frame
    for i in eachindex(points)
        v = points[i]-plane.point
        answer[i] = eltype(answer)(dot(o_x,v), dot(o_y,v), dot(o_z,v))
    end
    answer
end

"""
    compatiblesPlane(plane, points, normals, eps, alpharad)

Create a bool-indexer array for those points that are compatible to the plane.
Give back the projected points too for parameter space magic.

Compatibility is measured with an `eps` distance to the plane and an `alpharad` angle to it's normal.
"""
function compatiblesPlane(plane, points, normals, eps, alpharad)
    @assert length(points) == length(normals) "Size must be the same."
    projecteds = project2plane(plane, points)
    # eps check
    c1 = [abs(a[3]) < eps for a in projecteds]
    # alpha check
    c2 = [isparallel(plane.normal, normals[i], alpharad) && c1[i] for i in eachindex(normals)]
    # projecteds[c2] are the compatible points
    return c2, projecteds
end

"""
    bitmapparameters(parameters, compatibility, beta, idsource)

Bitmap the compatible parameters. An `idsource` is used to create the indexer map.

# Arguments:
- `parameters::AbstractArray`: contains the parameter pairs.
- `compatibility::AbstractArray`: indicates if the i-th parameter-pair is compatible.
- `beta<:Real`: resolution of the bitmap.
- `idsource::AbstractArray`: source of the indexer map.
"""
function bitmapparameters(parameters, compatibility, beta, idsource)
    @assert length(parameters) == length(compatibility) == length(idsource) "Everything must have the same length."
    miv, mav = findAABB(parameters)
    # perturbe minv and maxv
    pert = convert(typeof(miv), fill(0.1, length(miv)))
    minv = miv-pert
    maxv = mav+pert
    xs = round(Int, (maxv[1]-minv[1])/beta)
    ys = round(Int, (maxv[2]-minv[2])/beta)
    # TODO: ransac(pcr, αα, ϵϵ, tt, usegloββ, connekeyy, ptt, ττ, itermax)
    @assert xs > 0 && ys > 0 "max-min should be positive. Check the code! xs: $xs, ys:$ys"
    βx = (maxv[1]-minv[1])/xs
    βy = (maxv[2]-minv[2])/ys

    bitmap = falses(xs,ys)
    idxmap = zeros(Int,xs,ys)
    for i in eachindex(parameters)
        if compatibility[i]
            xplace = ceil(Int,(parameters[i][1]-minv[1])/βx)
            yplace = ceil(Int,(parameters[i][2]-minv[2])/βy)
            if xplace!=0 && yplace!=0 && xplace!=xs && yplace!=ys
                if !bitmap[xplace, yplace]
                    bitmap[xplace, yplace] = true
                    idxmap[xplace, yplace] = idsource[i]
                else
                    @warn "Project 2 or more into one! iteration:$i, id:$(idsource[i])."
                end
            else
                @warn "boundserror at $i: ($xplace,$yplace)!"
            end
        end
    end
    return bitmap, idxmap, (βx, βy)
end

"""
    bitmapparameters(parameters, compatibility, beta)

Bitmap the compatible parameters. The `idsource` is: `1:length(parameters)`.

# Arguments:
- `parameters::AbstractArray`: contains the parameter pairs.
- `compatibility::AbstractArray`: indicates if the i-th parameter-pair is compatible.
- `beta<:Real`: resolution of the bitmap.
"""
function bitmapparameters(parameters, compatibility, beta)
    return bitmapparameters(parameters, compatibility, beta, 1:length(parameters))
end

"""
    compatiblesSphere(plane, points, normals, eps, alpharad)

Create a bool-indexer array for those points that are compatible to the sphere.
Give back the projected points too for parameter space magic.
Return a bool indexer for (under,over) too.

Compatibility is measured with an `eps` distance to the sphere and an `alpharad` angle to it's normal.
"""
function compatiblesSphere(sphere, points, normals, eps, alpharad)
    @assert length(points) == length(normals) "Size must be the same."
    # eps check
    o = sphere.center
    R = sphere.radius
    c1 = [abs(norm(a-o)-R) < eps for a in points]
    # alpha check
    α = cos(alpharad)
    if sphere.outwards
        c2 = [isparallel(normalize(points[i]-o), normals[i], α) && c1[i] for i in eachindex(points)]
    else
        c2 = [isparallel(normalize(o-points[i]), normals[i], α) && c1[i] for i in eachindex(points)]
    end

    under = falses(length(points))
    over = falses(length(points))
    for i in eachindex(points)
        if points[i][3] <= o[3]
            under[i] = true
        else
            over[i] = true
        end
    end
    # TODO
    # normalize by it's own length or maximum length??? r+3ϵ
    proj_points = [SVector{2}(normalize(a[1:2]-o[1:2])) for a in points]
    param_points = unitdisk2square.(proj_points)
    return c2, (under=under, over=over), param_points
end

"""
    compatiblesCylinder(cylinder, points, normals, eps, alpharad)

Create a bool-indexer array for those points that are compatible to the cylinder.
Give back the projected points too for parameter space magic.

Compatibility is measured with an `eps` distance to the cylinder and an `alpharad` angle to it's normal.
"""
function compatiblesCylinder(cylinder, points, normals, eps, alpharad)
    @assert length(points) == length(normals) "Size must be the same."

    c = cylinder.center
    R = cylinder.radius
    a = cylinder.axis

    comp = falses(length(points))
    pars = fill(SVector(0.0,0.0), length(points))

    for i in 1:length(points)
        curr_norm = points[i] - a*dot( a, points[i]-c ) - c
        # if the radius is correct
        if abs(norm(curr_norm)-R) < eps
            if cylinder.outwards
                comp[i] = isparallel(normalize(curr_norm), normals[i], alpharad)
            else
                comp[i] = isparallel(-normalize(curr_norm), normals[i], alpharad)
            end
        end

        # playing with parameters:
        #TODO: implement it
    end
    return comp, pars
end

"""
    largestconncomp(bimage, indmap, connectivity)

Search for the largest connected component on a binary image.

Return the indexes from `indmap` that are part of the largest component.
"""
function largestconncomp(bimage, indmap, connectivity::AbstractArray)
    @assert size(bimage) == size(indmap) "Binary image and index map are not the same size!"
    # labeling
    labeled = label_components(bimage, connectivity)
    # size of each labelled area
    lab_length = component_lengths(labeled)
    # largest connected component (lab_length[1] is the background)
    # TODO: ERROR: ArgumentError: collection must be non-empty
    max_l = argmax(lab_length[2:end])+1
    # indexes that counts towards the largest area
    inds = findall(x->x==max_l-1, labeled)
    # get all the indexes that count towards the largest component
    return indmap[inds]
end

"""
    largestconncomp(bimage, indmap, conn_key = :default)

Search for the largest connected component on a binary image.

Return the indexes from `indmap` that are part of the largest component.

# Arguments:
- `bimage::BitArray{2}`: binary image, where `true` is the object.
- `indmap::AbstractArray`: index map, with the size of `bimage`.
- `conn_key::Symbol`: currently `:default` or `:eight` connectivity.
"""
function largestconncomp(bimage, indmap, conn_key::Symbol = :default)
    if conn_key == :default
        return largestconncomp(bimage, indmap, 1:ndims(bimage))
    elseif conn_key == :eight
        return largestconncomp(bimage, indmap, trues(3,3))
    else
        @error("No such key implemented: $conn_key.")
    end
end

"""
    refit(s, pc, ϵ, α)

Refit plane. Only s.inpoints is updated.
"""
function refitplane(s, pc, ϵ, α)
    # TODO: use octree for that
    cp, _ = compatiblesPlane(s.candidate.shape, pc.vertices[pc.isenabled], pc.normals[pc.isenabled], ϵ, α)
    s.inpoints = ((1:pc.size)[pc.isenabled])[cp]
    s
end

"""
    refit(s, pc, ϵ, α)

Refit sphere. Only s.inpoints is updated.
"""
function refitsphere(s, pc, ϵ, α)
    # TODO: use octree for that
    cpl, uo, sp = compatiblesSphere(s.candidate.shape, pc.vertices[pc.isenabled], pc.normals[pc.isenabled], ϵ, α)
    # verti: összes pont indexe, ami enabled és kompatibilis
    verti = (1:pc.size)[pc.isenabled]
    underEn = uo.under .& cpl
    overEn = uo.over .& cpl
    s.inpoints = count(underEn) >= count(overEn) ? verti[underEn] : verti[overEn]
    s
end

"""
    refit(s, pc, ϵ, α)

Refit cylinder. Only s.inpoints is updated.
"""
function refitcylinder(s, pc, ϵ, α)
    # TODO: use octree for that
    cp, _ = compatiblesCylinder(s.candidate.shape, pc.vertices[pc.isenabled], pc.normals[pc.isenabled], ϵ, α)
    s.inpoints = ((1:pc.size)[pc.isenabled])[cp]
    s
end
