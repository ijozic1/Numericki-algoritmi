using Pkg;
Pkg.add("NLsolve")

using NLsolve, LinearAlgebra,Plots;


function f!(F, x)
    F[1] = cos(x[1])-x[1]
end
nlsolve(f!, [2.1])
x = [-2:0.2:2;]
y=cos.(x);
plot(x,y);
plot!(x,x)

function f!(F, x)
    F[1] = x[2]-x[1] .^2+1
    F[2] = x[1] -(2*x[2]-x[2].^2)/3
end
x0 = [1.0,  0.0];
x1 = [-3.0 8.0];
nlsolve(f!, x0)
nlsolve(f!, x1)
x = [-3:0.2:3;];
y = x.^2 .+1;
plot(x,y);
y = [-3:0.2:5;];
x=(2 .*y-y.^2) ./3;
plot!(x,y)

Pkg.add("Optim")
using Optim;

f(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
optimize(f, [0.0, 0.0])

function rosenbrock(x::Vector)
    return (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
end
  
x, y = -1.5:0.1:1.5, -1.5:0.1:1.5
z = Surface((x,y)->rosenbrock([x,y]), x, y)
surface(x,y,z)


x = [ 1 3 5 4 2 ];
a = [ 1 5 4 2];
b = [ 4 1 3 ];
a1 = [ a 0 0 0 0 ];
b1 = [ b 0 0 0 0 0 ];
Pkg.add("FFTW"); Pkg.add("DSP")

using FFTW, DSP;
y=fft(x)
ifft(y)
y1=dct(x)
idct(y1)
conv(a,b)


Pkg.add("Polynomials"); using Polynomials;
p1=Polynomial([1,5,4,2])
p2=Polynomial([4,1,3])
p1*p2

a1f=fft(a1)
b1f=fft(b1)
c1f=a1f.*b1f
ifft(c1f)

#ZSR
Pkg.add("DifferentialEquations")
using DifferentialEquations
function opruga(du,u,p,t)
    du[1] = u[2];
    du[2] = -p[1]*u[1]/p[2]
end

p = [100.0;10.0]
u0 = [0.0;1.0]
tspan = (0.0,10.0)
prob = ODEProblem(opruga,u0,tspan,p)
sol = solve(prob)

plot(sol,vars=(1))
plot(sol,vars=(2))