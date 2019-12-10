# module ProFit

## References ##

# written by: Péter Salvi
# original source: https://github.com/salvipeter/profile-fitter

# The general workflow is based on:
#
# P. Benkő, T. Várady:
#   Best fit translational and rotational surfaces for reverse engineering shapes.
#   In: The mathematics of surfaces IX, pp. 70-81, 2000.

# The thinning algorithm for the guiding polygon is based on:
#
# I-K. Lee:
#   Curve reconstruction from unorganized points.
#   In: Computer Aided Geometric Design, Vol. 17(2), pp. 161-177, 2000.

# Circle fitting is based on:
#
# S.J. Ahn, W. Rauh, H-J. Warnecke:
#   Least-squares orthogonal distances fitting of circle, sphere, ellipse, hyperbola, and parabola.
#   In: Pattern Recognition, Vol. 34(12), pp. 2283-2303, 2001.

## TODO ##

# Beautification:
# - delete short arcs/segments
# - unite similar arcs
# - create segments from arcs with large radii
# - ensure C0 continuity
# - ensure G1 continuity where plausible

#using LinearAlgebra
#using Plots
#using Random

# Delaunay triangulation

"""
    bounding_triangle(points)

Returns a triangle that contains all points.
"""
function bounding_triangle(points)
    minx, maxx = extrema((x->x[1]).(points))
    miny, maxy = extrema((x->x[2]).(points))
    center = [(minx + maxx) / 2, (miny + maxy) / 2]
    r = norm(center - [minx, miny])
    [center - [sqrt(3) * r, r],
     center + [0, 2r],
     center + [sqrt(3) * r, -r]]
end

"""
    circumcircle(a, b, c)

Given three points, gives back the circumcircle as the tuple (center, radius).
"""
function circumcircle(a, b, c)
    d = 2 * (a[1] * (b[2] - c[2]) + b[1] * (c[2] - a[2]) + c[1] * (a[2] - b[2]))
    if abs(d) < 1e-10
        # Degenerate case
        if norm(a - c) > norm(a - b)
            b, c = c, b
        end
        return ((a + b) / 2, norm(a - b) / 2)
    end
    x = (a[1]^2 + a[2]^2) * (b[2] - c[2]) +
        (b[1]^2 + b[2]^2) * (c[2] - a[2]) +
        (c[1]^2 + c[2]^2) * (a[2] - b[2])
    y = (a[1]^2 + a[2]^2) * (c[1] - b[1]) +
        (b[1]^2 + b[2]^2) * (a[1] - c[1]) +
        (c[1]^2 + c[2]^2) * (b[1] - a[1])
    center = [x / d, y / d]
    center, norm(center - a)
end

"""
    to_edges!(triangles)

Converts a list of triangles (index-triples) to a list of edges (index-pairs).
The indices in the individual edges are sorted.
"""
function to_edges!(triangles)
    edges(tri) = sort.([[tri[1], tri[2]], [tri[1], tri[3]], [tri[2], tri[3]]])
    mapreduce(edges, append!, triangles)
end

"""
    walk_polygon(edges)

Given a list of edges (pairs of indices) defining a closed loop,
returns a sorted list of points going around the edges.
"""
function walk_polygon(edges)
    poly = [edges[1][1], edges[1][2]]
    while length(poly) != length(edges)
        for edge in edges
            if edge[1] == poly[end] && edge[2] != poly[end-1]
                push!(poly, edge[2])
                break
            elseif edge[2] == poly[end] && edge[1] != poly[end-1]
                push!(poly, edge[1])
                break
            end
        end
    end
    poly
end

"""
    delaunay(points)

Returns a list of triangles, each consisting of 3 indices to the input list.
"""
function delaunay(points)
    n = length(points)
    p = copy(points)
    append!(p, bounding_triangle(points))
    triangles = [[n+1, n+2, n+3]]

    for i in 1:n
        bad = filter(triangles) do tri
            c, r = circumcircle(p[tri[1]], p[tri[2]], p[tri[3]])
            norm(p[i] - c) <= r
        end
        setdiff!(triangles, bad)
        #isempty(bad) && continue
        bad_edges = to_edges!(bad)
        bad_boundary = filter(bad_edges) do edge
            count(x -> x == edge, bad_edges) == 1
        end
        poly = walk_polygon(bad_boundary)
        m = length(poly)
        for j in 1:m
            push!(triangles, [i, poly[j], poly[mod1(j+1,m)]])
        end
    end

    filter(triangles) do tri
        all(x -> x <= n, tri)
    end
end

# Thinning

"""
    spanning_tree(edges, weights)

Given a list of edges (index pairs) and associated weights,
returns a list of edges defining a minimal weight spanning tree.
"""
function spanning_tree(edges, weights)
    order = sortperm(weights)
    result = Array{Array{Int64,1},1}(undef, 0)
    sets = Dict()
    for edge in edges
        for p in edge
            sets[p] = p
        end
    end
    for i in order
        a, b = edges[i]
        if sets[a] != sets[b]
            push!(result, edges[i])
            new, old = minmax(sets[a], sets[b])
            for (k, v) in sets
                if v == old
                    sets[k] = new
                end
            end
        end
    end
    result
end

"""
    thinning(points, edges, thickness)

Given an Euclidean minimum spanning tree over a set of points taken from a (noisy) curve,
this function generates an ordered list of points representing a thinned version of the curve.
The `thickness` parameter should be around the same size as the noise.
"""
function thinning(points, edges, thickness)
    lookup = [Int[] for _ in points]
    for e in edges
        push!(lookup[e[1]], e[2])
        push!(lookup[e[2]], e[1])
    end

    function bfs(start)
        seen = fill(false, length(points))
        result = []
        queue = [start]
        while !isempty(queue)
            v = popfirst!(queue)
            push!(result, v)
            for w in lookup[v]
                if !seen[w]
                    push!(queue, w)
                    seen[w] = true
                end
            end
        end
        result
    end
    function neighbors(start)
        x = points[start]
        result = []
        function rec(next)
            push!(result, next)
            for p in lookup[next]
                if !(p in result) && norm(x - points[p]) < thickness
                    rec(p)
                end
            end
        end
        rec(start)
        result
    end

    ordered = bfs(bfs(1)[end])
    thinned = Array{Array{Float64,1},1}(undef, 0)
    chunks = Array{Array{Int64,1},1}(undef, 0)
    while !isempty(ordered)
        adjacent = neighbors(ordered[1])
        push!(thinned, mapreduce(x -> points[x], +, adjacent) / length(adjacent))
        push!(chunks, adjacent)
        setdiff!(ordered, adjacent)
    end
    thinned, chunks
end

"""
    thinning(points, thickness)

Given a set of unsorted points taken from a (noisy) curve,
this function generates an ordered list of points representing a thinned version of the curve.
The `thickness` parameter should be around the same size as the noise.
"""
function thinning(points, thickness)
    tris = delaunay(points)
    all_edges = to_edges!(tris)
    weights = [norm(points[e[1]] - points[e[2]]) for e in all_edges]
    tree = spanning_tree(all_edges, weights)
    thinning(points, tree, thickness)
end

function alle(l)
    ae = Vector{Vector{Int}}(undef, 0)
    for i in 1:l
        for j in i+1:l
            #i == j && continue
            push!(ae, [i,j])
        end
    end
    ae
end

"""
    thinning_slow(points, thickness)

Given a set of unsorted points taken from a (noisy) curve,
this function generates an ordered list of points representing a thinned version of the curve.
The `thickness` parameter should be around the same size as the noise.
Delaunay triangulation is not used to accelerate the code.
"""
function thinning_slow(points, thickness)
    @logmsg IterLow1 "Collect edges for $(length(points)) points."
    all_edges = alle(size(points,1))
    @logmsg IterLow1 "Computing distance"
    weights = [norm(points[e[1]] - points[e[2]]) for e in all_edges]
    @logmsg IterLow1 "Spanning tree"
    tree = spanning_tree(all_edges, weights)
    @logmsg IterLow1 "Thinning"
    thinning(points, tree, thickness)
end

function thinning_deldir(points, thickness)
    miv, mav = findAABB(points)
    rect = [miv[1]; mav[1]; miv[2]; mav[2]]
    pix = (x->x[1]).(points)
    piy = (x->x[2]).(points)
    @logmsg IterLow1  "Delaunay"
    deli1, deli2 = deldir(pix, piy, rect)
    all_edges = [ [deli1[i], deli2[i]] for i in 1:size(deli1,1)]
    #tris = delaunay(points)
    #all_edges = to_edges!(tris)
    @logmsg IterLow1 "Computing distance for $(length(all_edges)) edges, with $(length(points)) points"
    weights = [norm(points[e[1]] - points[e[2]]) for e in all_edges]
    @logmsg IterLow1 "Spanning tree"
    tree = spanning_tree(all_edges, weights)
    @logmsg IterLow1 "Thinning"
    thinning(points, tree, thickness)
end

# Fitting

"""
    fit_line(X)

Fits a least-squares line on the given points, and returns a `(point, normal)` pair,
where `normal` is a unit vector perpendicular to the line.
"""
function fit_line(X)
    centroid = sum(X) / length(X)
    A = reduce(hcat, [p - centroid for p in X])
    centroid, svd(A).U[:,2]
end

"""
    line_error(X, line)

Returns the squared sum of Euclidean distances to the line.
"""
function line_error(X, line)
    p, n = line
    sum(q -> dot(p - q, n)^2, X)
end

"""
    line_cut(X, line)

Cuts the segment relevant to the given points. Returns the endpoints.
"""
function line_cut(X, line)
    p, n = line
    project(q) = q + dot(p - q, n) * n
    project(X[1]), project(X[end])
end

"""
    fit_circle(X; tolerance, max_iteration, λ)

Fits a least-squares circle on the given points. `X` is a list of 2D points.
The iteration stops after the deviation between the last two approximation becomes
less than `tolerance`, or if the the iteration count exceeds `max_iteration`.
The step size for the Gauss-Newton iteration is given by `λ`.

The return value is a `(radius, center)` tuple.
"""
function fit_circle(X; tolerance = 1.0e-7, max_iteration = 1000, λ = 1.3)
    n = 2                       # dimension
    m = length(X)
    Xc = sum(X) / m
    R = sqrt(sum(Xi -> norm(Xi - Xc)^2, X) / m)
    a = vcat(R, Xc)
    J = Array{Float64}(undef, n * m, n + 1)
    rhs = Array{Float64}(undef, n * m)
    for _ in 1:max_iteration
        R = a[1]
        Xc = a[2:end]
        for i in 1:m
            v = X[i] - Xc
            d = norm(v)
            J[(i-1)*n+1:i*n,1] = v / d
            J[(i-1)*n+1:i*n,2:end] = I - R / d * (I - v * transpose(v) / d^2)
            rhs[(i-1)*n+1:i*n] = v / d * (d - R)
        end
        Δa = J \ rhs
        a += λ * Δa
        norm(Δa) < tolerance && return a[1], a[2:end]
    end
    # @warn "Circle fitting exited with maximum iteration"
    a[1], a[2:end]
end

"""
    circle_error(X, circle)

Returns the squared sum of Euclidean distances to the circle.
"""
function circle_error(X, circle)
    r, p = circle
    sum(q -> (norm(q - p) - r)^2, X)
end

"""
    circle_cut(X, circle)

Cuts the arc relevant to the given points. Returns the angles at the endpoints.
"""
function circle_cut(X, circle)
    p = circle[2]
    angle(q) = copysign(acos(normalize(q - p)[1]), asin(normalize(q - p)[2]))
    dist(a, b) = b > a ? b - a : 2pi + b - a
    from = angle(X[1])
    mid = angle(X[(length(X)+1)÷2])
    to = angle(X[end])
    if dist(from, mid) > dist(from, to)
        from, to = to, from
    end
    if to < from
        to += 2pi
    end
    from, to
end

"""
    fit_in_chunks(points, chunks, tolerance; seed_size)

Fits circular arcs and line segments on the given unordered points, with the help
of a guiding polygon.
The function also needs a list of `chunks` - an ordered collection of index lists
(the indices in each chunk are unordered, but the chunks progress through the curve).
Curve segments are separated when the fitting error is above tolerance.
All segment are assumed to consist of at least `seed_size` number of chunks.

The result is a list of named tuples of the following form:
  `(type=:arc, r=radius, p=center, min=minangle, max=maxangle)`
or
  `(type=:segment, a=start, b=end)`
"""
function fit_in_chunks(points, guide, chunks, tolerance; seed_size = 3)
    function bestfit(a, b)
        X = mapreduce(chunk -> points[chunk], append!, chunks[a:b])
        line = fit_line(X)
        circle = fit_circle(X)
        le, ce = line_error(X, line), circle_error(X, circle)
        if le < ce
            a, b = line_cut(guide[a:b], line)
            return le, (type=:segment, a=a, b=b)
        end
        a, b = circle_cut(guide[a:b], circle)
        ce, (type=:arc, r=circle[1], p=circle[2], min=a, max=b)
    end

    length(chunks) < seed_size && return [bestfit(1, length(chunks))[2]], [(1,length(chunks))]

    local best
    result = []
    result_chunks = []
    from = 1
    i = 1
    while i <= length(chunks) - seed_size + 1
        next = bestfit(i, i + seed_size - 1)[1]
        if from == i || next < best
            best = next
        elseif next - best > tolerance
            push!(result, bestfit(from, i + seed_size - 2)[2])
            push!(result_chunks, (from, i + seed_size - 2))
            from = i + seed_size - 2
            i = from - 1
        end
        i += 1
    end
    push!(result, bestfit(from, length(chunks))[2])
    push!(result_chunks, (from, length(chunks)))
    result, result_chunks
end

function fit(points, thickness)
    thinned, chunks = thinning(points, thickness)
    fit_in_chunks(points, thinned, chunks, thickness / 2)
end

function closeprofile!(segments)
    size(segments, 1) < 2 && return segments
    push!(segments, deepcopy(segments[1]))
    return segments
end

function centroid(points)
    com = sum(points)
    return com/size(points, 1)
end

# end # module
