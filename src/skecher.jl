using Pkg
Pkg.activate()

using Revise
using LinearAlgebra
using StaticArrays
using Makie

include("generatesamples.jl")
using .samples

# testing sphere parameter bitmap

using AbstractPlotting

o = SVector(0.0,5,π);
R = 13;
tsP, tsN = samplesphere(o, R, (4,4));
#sc =  Scene()
#showgeometry(sc, tsP, tsN)

parSp = Array{SArray{Tuple{2},Float64,1,2},1}()

for i in eachindex(tsP)
        vec = tsP[i]-o
        fí = atan(vec[2],vec[1])
        theta = atan(vec[3], vec[1])
        push!(parSp, SVector(fí, theta))
end

#sc2 = scatter(parSp)

hbox(scatter(tsP, markersize=2), scatter(parSp))

## Dating with Images
using Images
using ImageView

## basic test image
using TestImages
img = testimage("mandrill");
imshow(img);

## finding connected components on a binary image
# https://juliaimages.org/latest/function_reference/#Images.label_components
indexes = rand(Int,5,5);
tk = rand(Bool,5,5);
imshow(tk)

# labeling
lbs = label_components(tk);
# 8 connectivity
lbs8 = label_components(tk, trues(3,3))
# size of the components
lls = component_lengths(lbs)
# largest connected component (first is the background)
maxi = argmax(lls[2:end])+1
idxs = findall(x->x==maxi-1, lbs)
# get all the labels
lbs[idxs]
# get all the pixels
tk[idxs]
# get all the indexes
indexes[idxs]
# I don't know what is this for: idxs = component_indices(lbs)
