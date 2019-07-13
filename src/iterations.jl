function ransac(pc, αϵ, t, pt, τ, itmax, drawN, minleftover, setenabled; reset_rand = false)
    if setenabled
        pc.isenabled = trues(pc.size)
    end
    ransac(pc, αϵ, t, pt, τ, itmax, drawN, minleftover, reset_rand=reset_rand)
end

function ransac(pc, αϵ, t, pt, τ, itmax, drawN, minleftover; reset_rand = false)
    reset_rand && Random.seed!(1234)
    sp_ϵ = αϵ.ϵ_sphere
    sp_α = αϵ.α_sphere

    p_ϵ = αϵ.ϵ_plane
    p_α = αϵ.α_plane

    cy_ϵ = αϵ.ϵ_cylinder
    cy_α = αϵ.α_cylinder

    # build an octree
    @assert drawN > 2
    subsetN = length(pc.subsets)
    @debug "Building octree."
    minV, maxV = findAABB(pc.vertices);
    octree = Cell(SVector{3}(minV), SVector{3}(maxV), OctreeNode(pc, collect(1:pc.size), 1));
    r = OctreeRefinery(8);
    adaptivesampling!(octree, r);
    @debug "Octree finished."
    # initialize levelweight vector to 1/d
    # TODO: check if levelweight is not empty
    fill!(pc.levelweight, 1/length(pc.levelweight))
    fill!(pc.levelscore, zero(eltype(pc.levelscore)))

    random_points = randperm(pc.size);
    candidates = ShapeCandidate[]
    scoredshapes = ScoredShape[]
    extracted = ScoredShape[]
    # smallest distance in the pointcloud
    #lsd = smallestdistance(pc.vertices)
    @debug "Iteration begins."
    # iterate begin
    for k in 1:itmax
        if count(pc.isenabled) < minleftover
            @debug "Break at $k, because left only: $(length(findall(pc.isenabled)))"
            break
        end
        # generate t candidate
        for i in 1:t
            #TODO: that is unsafe, but probably won't interate over the whole pc
            # select a random point
            if length(random_points)<10
                random_points = randperm(pc.size)
                @debug "Recomputing randperm."
            end
            r1 = popfirst!(random_points)

            #TODO: helyettesíteni valami okosabbal,
            # pl mindig az első enabled - ha már nagyon sok ki van véve,
            # akkor az gyorsabb lesz
            while ! pc.isenabled[r1]
                r1 = rand(1:pc.size)
            end
            # search for r1 in octree
            current_leaf = findleaf(octree, pc.vertices[r1])
            # get all the parents
            cs = getcellandparents(current_leaf)
            # revese the order, cause it's easier to map with levelweight
            reverse!(cs)
            # chosse the level with the highest score
            # if multiple maximum, the first=largest cell will be used
            curr_level = argmax(pc.levelweight[1:length(cs)])
            # choose 3 more from cs[curr_level].data.incellpoints
            #sdf = shuffle(cs[curr_level].data.incellpoints)
            #an indexer array for random indexing
            cell_ind = cs[curr_level].data.incellpoints
            bool_cell_ind = pc.isenabled[cell_ind]
            # alábbi lehet view
            randomized_inds = cell_ind[bool_cell_ind]
            # nem shuffle hanem randperm kéne, mert itt már csak néhány random szám kell, mert mindegyik enabled
            shuffle!(randomized_inds)
            #println("cell_ind: $(size(cell_ind,1)), randomized_inds: $(size(randomized_inds,1))")
            sd = [r1]
            #kell egy ellenőrzés, hogy van-e drawN-1 db enabled pont
            #elég lehetne egy for ciklus
            for n in randomized_inds
                #ezen mindig végigiterálok, ki kell belőle szedni azt, ami nem enabled
                # don't use the same point twice
                n == r1 && continue
                # use only the enabled points
                #@show pc.isenabled[n]
                push!(sd, n)
                length(sd) == drawN && break
            end
            #=
            while length(sd) < drawN && length(sdf) > 0
                n = popfirst!(sdf)
                # don't use the same point twice
                n == r1 && continue
                # use only the enabled points
                pc.isenabled[n] && push!(sd, n)
            end
            =#
            # sd: 3 indexes of the actually selected points

            # if there's no 3 points, continue to the next draw
            length(sd) < drawN && continue

            #TODO: this should be something more general
            # fit plane to the selected points
            fp = isplane(pc.vertices[sd], pc.normals[sd], p_α)
            isshape(fp) && push!(candidates, ShapeCandidate(fp, curr_level))
            # fit sphere to the selected points
            sp = issphere(pc.vertices[sd], pc.normals[sd], sp_ϵ, sp_α)
            isshape(sp) && push!(candidates, ShapeCandidate(sp, curr_level))
            cp = iscylinder(pc.vertices[sd], pc.normals[sd], cy_ϵ, cy_α)
            isshape(cp) && push!(candidates, ShapeCandidate(cp, curr_level))
        end # for t

        # evaluate the compatible points, currently used as score
        # TODO: do something with octree levels and scores

        for c in candidates
            #TODO: save the bitmmaped parameters for debug
            which_ = 1
            ps = @view pc.vertices[pc.subsets[which_]]
            ns = @view pc.normals[pc.subsets[which_]]
            ens = @view pc.isenabled[pc.subsets[which_]]

            if isa(c.shape, FittedPlane)
                # plane
                # cp, pp = compatiblesPlane(c.shape, pc.vertices[pc.isenabled], pc.normals[pc.isenabled], ϵ, α)
                cp, pp = compatiblesPlane(c.shape, ps, ns, p_ϵ, p_α)
                inder = cp.&ens
                inpoints = (pc.subsets[which_])[inder]
                #inpoints = ((pc.subsets[1])[ens])[cp]
                score = estimatescore(length(pc.subsets[which_]), pc.size, length(inpoints))
                pc.levelscore[c.octree_lev] = pc.levelscore[c.octree_lev] + E(score)
                push!(scoredshapes, ScoredShape(c, score, inpoints))
            elseif isa(c.shape, FittedSphere)
                # sphere
                # cpl, uo, sp = compatiblesSphere(c.shape, pc.vertices[pc.isenabled], pc.normals[pc.isenabled], ϵ, α)
                cpl, uo, sp = compatiblesSphere(c.shape, ps, ns, sp_ϵ, sp_α)
                # verti: összes pont indexe, ami enabled és kompatibilis
                # lenne, ha működne, de inkább a boolean indexelést machináljuk
                verti = pc.subsets[1]
                underEn = uo.under .& cpl
                overEn = uo.over .& cpl

                inpoints = count(underEn) >= count(overEn) ? verti[underEn] : verti[overEn]
                score = estimatescore(length(pc.subsets[which_]), pc.size, length(inpoints))
                pc.levelscore[c.octree_lev] = pc.levelscore[c.octree_lev] + E(score)
                push!(scoredshapes, ScoredShape(c, score, inpoints))
            elseif isa(c.shape, FittedCylinder)
                cp, pp = compatiblesCylinder(c.shape, ps, ns, cy_ϵ, cy_α)
                inder = cp.&ens
                inpoints = (pc.subsets[which_])[inder]
                #inpoints = ((pc.subsets[1])[ens])[cp]
                score = estimatescore(length(pc.subsets[which_]), pc.size, length(inpoints))
                pc.levelscore[c.octree_lev] = pc.levelscore[c.octree_lev] + E(score)
                push!(scoredshapes, ScoredShape(c, score, inpoints))
            else
                # currently nothing else is implemented
                @warn "What the: $c"
            end # if
        end # for c
        # by now every candidate is scored into scoredshapes
        empty!(candidates)

        if length(scoredshapes) > 0
            # search for the largest score == length(inpoints) (for now)
            # best = largestshape(scoredshapes)
            best = findhighestscore(scoredshapes)
            bestshape = scoredshapes[best.index]
            # TODO: refine if best.overlap
            scr = E(bestshape.score)
            lengttt = length(bestshape.inpoints)
            ppp = prob(lengttt*subsetN, length(scoredshapes), pc.size, k=drawN)
            if k%50 == 0
                @debug "$k. it, best: $lengttt db, score: $scr, prob: $ppp, scored shapes: $(length(scoredshapes)) db."
            end
            #TODO: length will be only 1/numberofsubsets
            # if the probability is large enough, extract the shape
            if ppp > pt
                @debug "Extraction! best score: $(E(bestshape.score)), length: $(length(bestshape.inpoints))"

                # TODO: proper refit, not only getting the points that fit to that shape
                # what do you mean by refit?
                # refit on the whole pointcloud
                if bestshape.candidate.shape isa FittedPlane
                    refitplane(bestshape, pc, p_ϵ, p_α)
                elseif bestshape.candidate.shape isa FittedSphere
                    refitsphere(bestshape, pc, sp_ϵ, sp_α)
                elseif bestshape.candidate.shape isa FittedCylinder
                    refitcylinder(bestshape, pc, cy_ϵ, cy_α)
                else
                    @error "Whatt? panic with $(typeof(bestshape.candidate.shape))"
                end

                # invalidate points
                for a in bestshape.inpoints
                    pc.isenabled[a] = false
                end
                # extract the shape and delete from scoredshapes
                push!(extracted, deepcopy(bestshape))
                deleteat!(scoredshapes, best.index)
                # mark scoredshapes that have invalid points
                toremove = Int[]
                for i in eachindex(scoredshapes)
                    for a in scoredshapes[i].inpoints
                        if ! pc.isenabled[a]
                            push!(toremove, i)
                            break
                        end
                    end
                end

                # remove scoredshapes that have invalid points
                deleteat!(scoredshapes, toremove)
            end # if extract shape
        end # if length(scoredshapes)
        # update octree levels
        updatelevelweight(pc)

        # check exit condition
        # TODO: τ-t is le kéne osztani a subsestek számával
        if prob(τ/subsetN, length(scoredshapes), pc.size, k=drawN) > pt
            @debug "Break, at this point all shapes should be extracted: $k. iteráció."
            break
        end
        #if mod(k,itermax/10) == 0
        #    @debug "Iteration: $k"
        #end
    end # iterate end
    @debug "Iteration finished with $(length(extracted)) extracted and $(length(scoredshapes)) scored shapes."
    return scoredshapes, extracted
end # ransac function

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
