function testshape(shape, vertices, normals, epsilon, alpha)
    # compatiblesCylinder(cylinder, points, normals, params)
    rpars = ransacparameters([FittedCylinder], cylinder=(ϵ=epsilon, α=alpha,))
    compatiblesCylinder(shape, vertices, normals, rpars)
end
