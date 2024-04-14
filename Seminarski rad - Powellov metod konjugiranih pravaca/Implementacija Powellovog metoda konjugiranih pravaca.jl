using LinearAlgebra

function min_bracket(f, x0, hinit = 1e-5, hmax = 1e10, lambda = 1.4)
    h = hinit
    a = 0; b = 0
    f0 = 0; f1 = 0; f2 = 0

    if f(x0) > f(x0 + h) || f(x0) > f(x0 - h)
        while h < hmax
            a = x0 - h
            b = x0 + h
            f1 = f(a)
            f2 = f(b)
            f0 = f(x0)

            while isinf(f1)
                h /= 2 * (1 + lambda)
                a = x0 - h
                f1 = f(a)
            end

            while isinf(f2)
                h /= 2 * (1 + lambda)
                b = x0 + h
                f2 = f(b)
            end

            if f1 > f0 && f2 > f0
                return (a, b, x0)
            end

            h *= lambda

            if f0 > f1
                x0 = a
            elseif f0 > f2
                x0 = b
            end
        end

        throw(DomainError(-2, "Minimum could not be bracketed"))
    else
        return (x0 - h, x0 + h, x0)
    end
end

function golden_ratio_min(f, a, b, c, eps = 1e-8)
    phi = (1 + sqrt(5)) / 2
    d = 0

    if abs(c - a) < abs(b - c)
        d = b - (b - c) / phi
    else    
        d = c
        c = a + (c - a) / phi
    end

    u = f(c)
    v = f(d)

    while abs(b - a) > eps
        if u < v
            b = d
            d = c
            c = a + (c - a) / phi
            v = u
            u = f(c)
        else
            a = c
            c = d
            d = b - (b - d) / phi
            u = v
            v = f(d)
        end
    end

    return (a + b) / 2
end

function find_minimum(f)
    (a, b, c) = (0.0, 0.0, 0.0)
    try
        (a, b, c) = min_bracket(f, 0)
    catch
        throw(DomainError(-3, "Minimum could not be bracketed"))
    end
     
    return golden_ratio_min(f, a, b, c)
end

function powell_min(f, x0, points, maxiter = 100, eps = 1e-8)
    n = length(x0)
    x = zeros(n, n)
    u = I(n) + zeros(n, n)
    h = 0

    for k in 1 : maxiter
        push!(points, x0)

        try
            h = find_minimum((h) -> f(x0 + h * u[:, 1]))
        catch
            break
        end
        
        x[:, 1] = x0 + h * u[:, 1]

        for i in 2 : n
            try
                h = find_minimum((h) -> f(x[:, i - 1] + h * u[:, i]))
            catch
                break
            end
            
            x[:, i] = x[:, i - 1] + h * u[:, i]
        end

        for i in 1 : n - 1
            u[:, i] = u[:, i + 1]
        end

        u[:, n] = x[:, n] - x0
        if norm(u[:, n]) < eps
            return x0
        end

        try
            h = find_minimum((h) -> f(x[:, n] + h * u[:, n]))
        catch
            break
        end
        
        x0 = x[:, n] + h * u[:, n]
    end

    throw(DomainError(-1, "Minimum has not been found"))
end


#=testne funkcije=#
using PlotlyJS, Plots

plotlyjs()

#quadratic
function quadratic(x)
    return x[1]^2 + 5*x[2]^2 + x[1]*x[2] - x[1] + x[2]
end
x_q = [-1.0, -1.0]
points = []

try
    x = powell_min(quadratic, x_q, points)
    println(x)
    println(quadratic(x))
catch
    println("Exception")
end

X = []; Y = []
for i in eachindex(points)
    push!(X, points[i][1])
    push!(Y, points[i][2])
end

x, y = -3:0.1:3, -3:0.1:3
z = Surface((x,y) -> quadratic([x,y]), x, y)
fvector = quadratic.(points)

Plots.plot(X, Y, fvector, st=:scatter); 
plot!(x, y, z, st=:surface, c=:thermal)

points

#Rosenbrock
function rosenbrock(x)
    return (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
end

x_r2 = [2.0; 1.3]
points = []

try
    x = powell_min(rosenbrock, x_r2, points)
    println(x)
    println(rosenbrock(x))
catch
    println("Exception")
end

X = []; Y = []
for i in eachindex(points)
    push!(X, points[i][1])
    push!(Y, points[i][2])
end

x, y = -2.5:0.1:2.5, -2.5:0.1:2.5
z = Surface((x,y) -> rosenbrock([x,y]), x, y)
fvector = rosenbrock.(points)

Plots.plot(X, Y, fvector, st=:scatter); 
plot!(x, y, z, st=:surface, c=:thermal)

points

#Ackley
function ackley(x)
    term1 = -20 * exp(-0.2 * sqrt((x[1]^2 + x[2]^2) / 2))
    term2 = exp(0.5 * (cos(2 * x[1]) + cos(2 * x[2])))
    result = term1 - term2 + 20 + exp(1)
    return result
end
x_a = [1.0, 2.0]
points = []

try
    x = powell_min(ackley, x_a, points)
    println(x)
    println(ackley(x))
catch
    println("Exception")
end

X = []; Y = []
for i in eachindex(points)
    push!(X, points[i][1])
    push!(Y, points[i][2])
end

x, y = -3.5:0.1:3.5, -3.5:0.1:3.5
z = Surface((x,y) -> ackley([x,y]), x, y)
fvector = ackley.(points)

Plots.plot(X, Y, fvector, st=:scatter); 
plot!(x, y, z, st=:surface, c=:thermal)

points 

#Matyas
function matyas(x)
    return 0.26*(x[1]^2 + x[2]^2) - 0.48*x[1]*x[2]
end
x_m = [-1.0, -1.0]
points = []

try
    x = powell_min(matyas, x_m, points)
    println(x)
    println(matyas(x))
catch
    println("Exception")
end

X = []; Y = []
for i in eachindex(points)
    push!(X, points[i][1])
    push!(Y, points[i][2])
end

x, y = -2.5:0.1:2.5, -2.5:0.1:2.5
z = Surface((x,y) -> matyas([x,y]), x, y)
fvector = matyas.(points)

Plots.plot(X, Y, fvector, st=:scatter); 
plot!(x, y, z, st=:surface, c=:thermal)

points

#Himmelblau
function himmelblau(x)
    return (x[1]^2 + x[2]-11)^2 + (x[1] + x[2]^2 - 7)^2
end
x_h = [0.0, 0.0]
points = []

try
    x = powell_min(himmelblau, x_h, points)
    println(x)
    println(himmelblau(x))
catch
    println("Exception")
end

X = []; Y = []
for i in eachindex(points)
    push!(X, points[i][1])
    push!(Y, points[i][2])
end

x, y = -4.5:0.1:4.5, -4.5:0.1:4.5
z = Surface((x,y) -> himmelblau([x,y]), x, y)
fvector = himmelblau.(points)

Plots.plot(X, Y, fvector, st=:scatter); 
plot!(x, y, z, st=:surface, c=:thermal)

points

#Sfera
function sphere(x)
    return x[1]^2 + x[2]^2
end
x_s = [3.0, 1.7]
points = []

try
    x = powell_min(sphere, x_s, points)
    println(x)
    println(sphere(x))
catch
    println("Exception")
end

X = []; Y = []
for i in eachindex(points)
    push!(X, points[i][1])
    push!(Y, points[i][2])
end

x, y = -2.5:0.1:2.5, -2.5:0.1:2.5
z = Surface((x,y) -> sphere([x,y]), x, y)
fvector = sphere.(points)

Plots.plot(X, Y, fvector, st=:scatter); 
plot!(x, y, z, st=:surface, c=:thermal)

points

#proizvoljna f-ja vise od 2 varijable
function f5(x)
    return 2 * sin(x[1]^2 + 2*x[2]^2 + 10*x[3]^2 + x[4]^2 + 12*x[5]^2 - 5) + 3
end

x05 = [1.0; 2.0; 0.3; 3.3; 1.2]
points = []

try
    x = powell_min(f5, x05, points)
    println(x)
    println(f5(x))
catch
    println("Exception")
end

points 