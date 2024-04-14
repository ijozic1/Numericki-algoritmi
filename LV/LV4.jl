using Pkg; 
Pkg.add("Polynomials"); 
using Plots, Polynomials;

p1=Polynomial([1:2:5;])
p2=fromroots([1:3:10;])
roots(p1)
roots(p2)
coeffs(p1)
coeffs(p2)

p1+p2
p1-p2
p1*p2
p2 ÷ p1

p=Polynomial([1,2,3],:y)
p(5)
p(5+2*im)
p.([1:5;])

x=variable()
x^2

cheb=ChebyshevT([1, 0, 3, 4])
p = convert(Polynomial, cheb)
domain(p)
x = -10:10
extrema(mapdomain(ChebyshevT, x))

xs = range(0, 10, length=10)
ys = @. exp(-xs)
f = fit(xs, ys) # degree = length(xs) - 1
f2 = fit(xs, ys, 2) # degree = 2

scatter(xs, ys, markerstrokewidth=0, label="Data")
plot!(f, extrema(xs)..., label="Fit")
plot!(f2, extrema(xs)..., label="Quadratic Fit")

using Images,TestImages,ImageMagick

path = "coins.png"
img = load(path)
img_gray=Gray.(img)
img_mat=convert(Array{Float64}, img_gray)
edges, counts = imhist(img_mat,255)
histogram(vec(reinterpret(UInt8, img_mat)))
scatter( edges,counts[1:end-1], label="points")
f = fit(edges, counts[1:end-1],30)
y=f.(edges)
plot!(edges,y)
V=sort(abs.(diff(y)))
ind=sortperm(abs.(diff(y)))
secondd=diff(V)
valuem,i=findmin(secondd)
tresh=ind[i]
valueTresh = y[tresh]/255
img_gray
A = convert(Array{Float64}, img_gray);
img_binary = (A.> valueTresh).*A;
Gray.(img_binary)