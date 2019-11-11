"""
    bitmapparameters(parameters, compatibility, beta, idsource)

Bitmap the compatible parameters. An `idsource` is used to create the indexer map.

# Arguments:
- `parameters::AbstractArray`: contains the parameter pairs.
- `compatibility::AbstractArray`: indicates if the i-th parameter-pair is compatible.
- `beta<:Real`: resolution of the bitmap.
- `idsource::AbstractArray`: source of the indexer map.
"""
function bitmapparameters(parameters, compatibility, params, idsource)
    @unpack β = params
    @assert length(parameters) == length(compatibility) == length(idsource) "Everything must have the same length."
    miv, mav = findAABB(parameters)
    # perturbe minv and maxv
    pert = convert(typeof(miv), fill(0.1, length(miv)))
    minv = miv-pert
    maxv = mav+pert
    xs = round(Int, (maxv[1]-minv[1])/β)
    ys = round(Int, (maxv[2]-minv[2])/β)
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
function bitmapparameters(parameters, compatibility, params)
    return bitmapparameters(parameters, compatibility, params, 1:length(parameters))
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
    subscripts = component_subscripts(labeled)[max_l]
    inds = Int[]
    for i in 1:size(subscripts, 1)
        ii = subscripts[i]
        append!(inds, indmap[ii[1], ii[2]])
    end
    # get all the indexes that count towards the largest component
    return inds
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
function largestconncomp(bimage, indmap, conn_key::Symbol=:default)
    if conn_key == :default
        return largestconncomp(bimage, indmap, 1:ndims(bimage))
    elseif conn_key == :eight
        return largestconncomp(bimage, indmap, trues(3,3))
    else
        @error("No such key implemented: $conn_key.")
    end
end
