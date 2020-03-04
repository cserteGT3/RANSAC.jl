"""
    rodrigues(nv, Θ)

Create a rotation matrix from a normalized axis and an angle (in radian).
Near-zero elements will be chopped to zero.
"""
function rodrigues(nv, Θ)
    et = eltype(nv)
    R = zeros(et,3,3)
    R = nv*nv' + cos(Θ).*(Matrix{et}(I, 3,3) - nv*nv') + sin(Θ).*crossprodtensor(nv)
    return R
end

"""
    crossprodtensor(v)

Create the cross product tensor of 3 dimensional vector.
"""
function crossprodtensor(v)
    [0 -v[3] v[2];
    v[3] 0 -v[1];
    -v[2] v[1] 0]
end

"""
    rodriguesrad(nv, ϑ)

Create a rotation matrix from an axis (`nv`) and an angle (`ϑ`) in radians.
"""
function rodriguesrad(nv, ϑ)
    nvn = normalize(nv)
    return rodrigues(nvn, ϑ)
end

"""
    rodriguesdeg(nv, ϑ)

Create a rotation matrix from an axis (`nv`) and an angle (`ϑ`) in degrees.
"""
function rodriguesdeg(nv, ϑ)
    nvn = normalize(nv)
    return rodrigues(nvn, deg2rad(ϑ))
end

## Collection of useful functions.


"""
    arbitrary_orthogonal(vec)

Create an arbitrary orthogonal vector to `vec`.
"""
function arbitrary_orthogonal(vec)
    @assert size(vec,1) == 3 "Implemented only for 3 dimensional."
    v = normalize(vec)
    b0 = (v[1] <  v[2]) && (v[1] <  v[3]);
    b1 = (v[2] <= v[1]) && (v[2] <  v[3]);
    b2 = (v[3] <= v[1]) && (v[3] <= v[2]);
    rv = normalize!([convert(Int, b0), convert(Int, b1), convert(Int, b2)])
    return convert(typeof(vec), cross(v, rv))
end

"""
    arbitrary_orthogonal2(vec)

Create an arbitrary orthogonal vector to `vec`.
"""
function arbitrary_orthogonal2(vec)
    @assert size(vec,1) == 3 "Implemented only for 3 dimensional."
    v = normalize(vec)
    rv = normalize(rand(3))
    return convert(typeof(vec), cross(v, rv))
end


"""
    isparallel(v1, v2, alpharad)

Check if `v1` and `v2` are parallel in an `alpharad` angle and point towards the same direction.

The vectors considered to be normalized.
"""
function isparallel(v1, v2, alpharad)
    dot(v1, v2) > cos(alpharad)
end


"""
    findAABB(points)

Find the axis aligned minimum bounding box of the given pointcloud.
"""
function findAABB(points)
    minV = convert(Array, deepcopy(points[1]))
    maxV = convert(Array, deepcopy(points[1]))
    for i in 1:length(points)
        a = points[i]
        for j in 1:length(minV)
            minV[j] = minV[j] > a[j] ? a[j] : minV[j]
            maxV[j] = maxV[j] < a[j] ? a[j] : maxV[j]
        end
    end
    return convert(eltype(points), minV), convert(eltype(points), maxV)
end


"""
    chopzaxis(points)

Chop the 3rd element of every point.
"""
function chopzaxis(points)
    @assert length(points[1]) == 3 "Implemented only for 3 long vectors."
    [SVector(a[1], a[2]) for a in points]
end


"""
    unitdisk2square(p)

Transform point on unit disk to unit sqaure.

Maps a unit circle to a unit sqaure [0,1]^2.
Source: shirley1997 - A low distortion map between disk and square
"""
function unitdisk2square(p)
    r = norm(p)
    phi = atan(p[2], p[1])
    if phi < -π/4
        phi += 2*π # in range [-pi/4,7pi/4]
    end
    if phi < π/4 # region 1
        a = r
        b = phi*a/(π/4)
    elseif phi < 3*π/4 # region 2
        b = r
        a = -(phi-π/2)*b/(π/4)
    elseif phi < 5*π/4 # region 3
        a = -r
        b = (phi - π)*a/(π/4)
    else # region 4
        b = -r
        a = -(phi-3*π/2)*b/(π/4)
    end
    # TODO:
    # decision: squre: [-1,1]^2 or [0,1]^2
    convert(typeof(p), [(a+1)/2, (b+1)/2])
end

"""
    smallestdistance(points)

Find the smallest distance between the points.
"""
function smallestdistance(points)
    @assert length(points) > 1 "At least two point is needed for that."
    ld = norm(points[2]-points[1])
    for i in eachindex(points)
        for j in eachindex(points)
            if i!=j
                d = norm(points[i]-points[j])
                ld = d < ld ? d : ld
            end
        end
    end
    ld
end

"""
    minmaxdistance(points)

Find the smallest, largest and average distance between the points.
"""
function minmaxdistance(points)
    @assert length(points) > 1 "At least two point is needed for that."
    ld = norm(points[2]-points[1])
    maxid = norm(points[2]-points[1])
    sumd = 0.
    sumi = 0
    for i in eachindex(points)
        for j in eachindex(points)
            if i!=j
                d = norm(points[i]-points[j])
                sumd += d
                sumi += 1
                ld = d < ld ? d : ld
                maxid = d > maxid ? d : maxid
            end
        end
    end
    avg_d = sumd/sumi
    return (mind=ld, maxd=maxid, avgd=avg_d)
end

"""
    not0minmaxdistance(points)

Find the smallest and largest distance between the points.
Zero distance is not taken into account.
"""
function not0minmaxdistance(points)
    @assert length(points) > 1 "At least two point is needed for that."
    ld = norm(points[2]-points[1])
    maxid = norm(points[2]-points[1])
    for i in eachindex(points)
        for j in eachindex(points)
            if i!=j
                d = norm(points[i]-points[j])
                d == 0 && continue
                ld = d < ld ? d : ld
                maxid = d > maxid ? d : maxid
            end
        end
    end
    return (mind=ld, maxd=maxid)
end

"""
    prob(n, s, N, k)

The probability of successful detection of a shape sized `n`,
from a point cloud size of `N`, with `k` size of minimal sets and `s` draws.

# Arguments:
- `n`: size of the shape.
- `s`: number of candidates that have been drawn.
- `N`: size of the point cloud.
- `k`: size of the minimal set required to define a shape candidate.
"""
prob(n, s, N, k) = 1-(1-(n/N)^k)^s


"""
    havesameelement(A, B)

Return true if `A` and `B` have at least one common element.
"""
function havesameelement(A, B)
    for a in A
        for b in B
            a == b && return true
        end
    end
    false
end

@with_kw struct RANSACParameters{R<:Real} @deftype R
    ϵ_plane = 0.3
    α_plane = deg2rad(5)

    ϵ_sphere = 0.3
    α_sphere = deg2rad(5)

    ϵ_cylinder = 0.3
    α_cylinder = deg2rad(5)

    ϵ_cone = 0.3
    α_cone = deg2rad(5)
    minconeopang = deg2rad(2)

    ϵ_transl = 0.3
    α_transl = deg2rad(5)
    # max deviaton of the normal being perpendicular to the translation direction
    α_perpend = cosd(89)

    ϵ_torus = 0.3
    α_torus = deg2rad(5)

    # number of points to be sampled (length of a minimal subset)
    drawN::Int = 3; @assert drawN>2
    # number of minimal sets sampled in one iteration
    minsubsetN::Int = 15; @assert minsubsetN>0
    # probability of detection
    prob_det = 0.9
    # minimal shape size
    τ::Int = 900
    # maximum number of iteration
    itermax::Int = 1000
    # if the number of enabled points fall under `leftover`,
    # the iteration terminates
    leftovers::Int = 1

    # threshold of two vectors being parallel (in degrees)
    parallelthrdeg = 1
    # threshold of points being collinear
    collin_threshold = 0.2
    # bitmap resolution
    β = 1
    # parameter in sphere fitting
    sphere_par = 0.02

    ## translational fitting parameters
    # ???
    diagthr = 0.1
    # conntectivity on the bitmap - not used currently
    transl_conn::Symbol = :eight
    # maximum number of contours on a plane
    #TODO: delete this
    max_group_num::Int = 3
    # maximum number of iterations of tryíng
    # to find < max_group_num number for contour patches
    #TODO: delete this
    max_contour_it::Int = 5
    thinning_par = 2.0
    # minimum % of the normals must be the same
    min_normal_num = 0.9
    # extract translational surface even though normals are not ok
    force_transl::Bool = false
    # thinning method: :slow/:fast/:deldir
    thin_method::Symbol = :slow
    # how close must they be to consider them as the same point?
    samep = Float64(eps(Float32))
    # check side parameter
    checksidepar = 0.04
    # "disabled by default"
    max_end_d = 10000.0
    # skip or jumpback
    jumpback::Bool = false

    # shapes that are fitted to the point cloud
    #shape_types::Array{Symbol,1} = [:sphere, :plane, :cylinder, :cone, :translational_surface]
    shape_types::Array{Symbol,1} = [:sphere, :plane, :cylinder, :cone]

    # track the number of candidates for the probabilities
    # 1. :lengthC - number of candidates in a given iteration
    # 2. :allcand - number of candidates that have been ever scored
    # 3. :nofminset - number of minimal sets that have been drawn so far
    # which to use for extracting the best candidate
    extract_s::Symbol = :nofminset
    # which to use to terminate the algorithm
    terminate_s::Symbol = :nofminset
end

"""
    allisdifferent(a::AbstractArray{T, 1}) where {T}

Check if all the elements of an array are different.
Supports only one dimensional arrays.
"""
function allisdifferent(a::AbstractArray{T, 1}) where {T}
    #TODO: make it more abstract with:
    #https://docs.julialang.org/en/v1/base/iterators/
    size(a, 1) < 2 && return true
    for i in 2:size(a, 1)
        for j in 1:i-1
            a[i] == a[j] && return false
        end
    end
    return true
end

"""
    setepsilons(p::RANSACParameters, ϵ)

Set all `ϵ_*` fields to `ϵ`.
"""
function setepsilons(p::RANSACParameters, ϵ)
    RANSACParameters(p, ϵ_plane=ϵ, ϵ_sphere=ϵ, ϵ_cylinder=ϵ, ϵ_cone=ϵ)
end

"""
    setalphas(p::RANSACParameters, α)

Set all `α_*` fields to `α`
"""
function setalphas(p::RANSACParameters, α)
    RANSACParameters(p, α_plane=α, α_sphere=α, α_cylinder=α, α_cone=α)
end

function chooseS(A, k)
    symtoint = Dict(:lengthC=>1, :allcand=>2, :nofminset=>3)
    return A[symtoint[k]]
end

#=
"""
    p2table(p)

Format a parameter into latex table.
"""
function p2table(p, fname, pre=true)
    b = IOBuffer()
    if pre
        println(b, """
        \\begin{table}[h]
        \\centering
        \\caption{Parameters for the segmentation for \\emph{Example ?}. (description of values: \\ref{tab:ransacpars})}
        """)
    end
    println(b,"""	\\begin{tabular}{ c | c | c | c }
		Parameter & Description & Parameter & Description \\\\
		\\hline
		\\hline
		\$N\$ & \$?\$ & \$p_t\$ & \$$(p.prob_det)\%\$\\\\
		\$\\tau\$ & \$$(p.τ)\$ & \$k\$ & \$3\$ \\\\
		\$d\$ & \$?\$& \$t\$ & \$$(p.minsubsetN)\$ \\\\
		\$r\$ & \$?\$ & - & - \\\\
		\$\\epsilon_{plane}\$ & \$$(p.ϵ_plane)\$ & \$\\alpha_{plane}\$  & \$$(p.α_plane)^{\\circ}\$\\\\
		\$\\epsilon_{sphere}\$ & \$$(p.ϵ_sphere)\$ & \$\\alpha_{sphere}\$  & \$$(p.α_sphere)^{\\circ}\$\\\\
		\$\\epsilon_{cylinder}\$ & \$$(p.ϵ_cylinder)\$ & \$\\alpha_{cylinder}\$  & \$$(p.α_cylinder)^{\\circ}\$\\\\
		\$\\epsilon_{cone}\$ & \$$(p.ϵ_cone)\$ & \$\\alpha_{cone}\$  & \$$(p.α_cone)^{\\circ}\$\\\\
        \$\\epsilon_{translational}\$ & \$$(p.ϵ_transl)\$ & \$\\alpha_{translational}\$  & \$$(p.α_transl)^{\\circ}\$\\\\
        \\end{tabular}""")
    if pre
        println(b, """	\\label{tab:ex1pars}
                    \\end{table}""")
    end
    open(fname, "w") do io
        println(io, b)
    end
end
=#
