using Pkg
Pkg.add(url="https://github.com/PumasAI/DataInterpolations.jl")
Pkg.add("Interpolations")
Pkg.add("Plots")

using DataInterpolations, Interpolations, Plots;

using LinearAlgebra;

x1=[LinRange(0,2*pi,100);]
y1=sin.(x1)
x2=[LinRange(-1,1,100);]
y2=1 ./ (1 .+ 25*x2.^2)
plot1=plot(x1,y1,title="Sinus");
plot2=plot(x2,y2,title="Runge");
plot(plot1,plot2,)

xdata1=[LinRange(0,2*pi,10);]
ydata1=sin.(xdata1);
xdata2=[LinRange(-1,1,10);]
ydata2=1 ./ (1 .+ 25*xdata2.^2);

A = DataInterpolations.LinearInterpolation(ydata1,xdata1)
scatter(xdata1,ydata1, label="Početne tačke")
xx1=[LinRange(0,2*pi,50);]
yy1=A.(xx1)
plot!(xx1,yy1,label="Interpolirane")
plot!(x1,y1,label="Stvarne")

B = DataInterpolations.LinearInterpolation(ydata2,xdata2)
scatter(xdata2,ydata2, label="Početne tačke")
xx2=[LinRange(-1,1,50);]
yy2=B.(xx2)
plot!(xx2,yy2,label="Interpolirane")
plot!(x2,y2,label="Stvarne")

A = DataInterpolations.QuadraticInterpolation(ydata1,xdata1)
scatter(xdata1,ydata1, label="Početne tačke")
xx1=[LinRange(0,2*pi,50);]
yy1=A.(xx1)
plot!(xx1,yy1,label="Interpolirane")
plot!(x1,y1,label="Stvarne")

B = DataInterpolations.QuadraticInterpolation(ydata2,xdata2)
scatter(xdata2,ydata2, label="Početne tačke")
xx2=[LinRange(-1,1,50);]
yy2=B.(xx2)
plot!(xx2,yy2,label="Interpolirane")
plot!(x2,y2,label="Stvarne")


A = DataInterpolations.ConstantInterpolation(ydata1,xdata1)
scatter(xdata1,ydata1, label="Početne tačke")
xx1=[LinRange(0,2*pi,50);]
yy1=A.(xx1)
plot!(xx1,yy1,label="Interpolirane")
plot!(x1,y1,label="Stvarne")

B = DataInterpolations.ConstantInterpolation(ydata2,xdata2)
scatter(xdata2,ydata2, label="Početne tačke")
xx2=[LinRange(-1,1,50);]
yy2=B.(xx2)
plot!(xx2,yy2,label="Interpolirane")
plot!(x2,y2,label="Stvarne")


A = DataInterpolations.LagrangeInterpolation(ydata1,xdata1)
scatter(xdata1,ydata1, label="Početne tačke")
xx1=[LinRange(0,2*pi,50);]
yy1=A.(xx1)
plot!(xx1,yy1,label="Interpolirane")
plot!(x1,y1,label="Stvarne")

B = DataInterpolations.LagrangeInterpolation(ydata2,xdata2)
scatter(xdata2,ydata2, label="Početne tačke")
xx2=[LinRange(-1,1,50);]
yy2=B.(xx2)
plot!(xx2,yy2,label="Interpolirane")
plot!(x2,y2,label="Stvarne")

A = DataInterpolations.CubicSpline(ydata1,xdata1)
scatter(xdata1,ydata1, label="Početne tačke")
xx1=[LinRange(0,2*pi,50);]
yy1=A.(xx1)
plot!(xx1,yy1,label="Interpolirane")
plot!(x1,y1,label="Stvarne")

B = DataInterpolations.CubicSpline(ydata2,xdata2)
scatter(xdata2,ydata2, label="Početne tačke")
xx2=[LinRange(-1,1,50);]
yy2=B.(xx2)
plot!(xx2,yy2,label="Interpolirane")
plot!(x2,y2,label="Stvarne")

A = DataInterpolations.BSplineInterpolation(ydata1,xdata1,2,:Uniform,:Uniform)
scatter(xdata1,ydata1, label="Početne tačke")
xx1=[LinRange(0,2*pi,50);]
yy1=A.(xx1)
plot!(xx1,yy1,label="Interpolirane")
plot!(x1,y1,label="Stvarne")

B = DataInterpolations.BSplineInterpolation(ydata2,xdata2,2,:Uniform,:Uniform)
scatter(xdata2,ydata2, label="Početne tačke")
xx2=[LinRange(-1,1,50);]
yy2=B.(xx2)
plot!(xx2,yy2,label="Interpolirane")
plot!(x2,y2,label="Stvarne")

itp = scale(interpolate(ydata1,BSpline(Linear()), OnGrid()),xdata1[1]:(xdata1[2]-xdata1[1]):xdata1[end])
xx1=[LinRange(-1pi,3pi,100);];
etpf = extrapolate(itp,Throw())
yy1=etpf(xx1)
scatter(xdata1,ydata1, label="Početne tačke")
plot!(xx1,yy1,label="Ekstrapolirane")
plot!(x1,y1,label="Stvarne")

itp = scale(interpolate(ydata2,BSpline(Linear()), OnGrid()),LinRange(-1,1,10))
xx2=[LinRange(-3,3,100);];
etpf = extrapolate(itp,Throw())
yy2=etpf(xx2)
scatter(xdata2,ydata2, label="Početne tačke")
plot!(xx2,yy2,label="Ekstrapolirane")
plot!(x2,y2,label="Stvarne")

itp = scale(interpolate(ydata1,BSpline(Linear()), OnGrid()),xdata1[1]:(xdata1[2]-xdata1[1]):xdata1[end])
xx1=[LinRange(-1pi,3pi,100);];
etpf = extrapolate(itp,Interpolations.Flat())
yy1=etpf(xx1)
scatter(xdata1,ydata1, label="Početne tačke")
plot!(xx1,yy1,label="Ekstrapolirane")
plot!(x1,y1,label="Stvarne")


itp = scale(interpolate(ydata2,BSpline(Linear()), OnGrid()),LinRange(-1,1,10))
xx2=[LinRange(-3,3,100);];
etpf = extrapolate(itp,Interpolations.Flat())
yy2=etpf(xx2)
scatter(xdata2,ydata2, label="Početne tačke")
plot!(xx2,yy2,label="Ekstrapolirane")
plot!(x2,y2,label="Stvarne")

itp = scale(interpolate(ydata1,BSpline(Linear()), OnGrid()),xdata1[1]:(xdata1[2]-xdata1[1]):xdata1[end])
xx1=[LinRange(-1pi,3pi,100);];
etpf = extrapolate(itp,Interpolations.Line())
yy1=etpf(xx1)
scatter(xdata1,ydata1, label="Početne tačke")
plot!(xx1,yy1,label="Ekstrapolirane")
plot!(x1,y1,label="Stvarne")

itp = scale(interpolate(ydata2,BSpline(Linear()), OnGrid()),LinRange(-1,1,10))
xx2=[LinRange(-3,3,100);];
etpf = extrapolate(itp,Interpolations.Line())
yy2=etpf(xx2)
scatter(xdata2,ydata2, label="Početne tačke")
plot!(xx2,yy2,label="Ekstrapolirane")
plot!(x2,y2,label="Stvarne")

itp = scale(interpolate(ydata1,BSpline(Linear()), OnGrid()),xdata1[1]:(xdata1[2]-xdata1[1]):xdata1[end])
xx1=[LinRange(-1pi,3pi,100);];
etpf = extrapolate(itp,Interpolations.Periodic())
yy1=etpf(xx1)
scatter(xdata1,ydata1, label="Početne tačke")
plot!(xx1,yy1,label="Ekstrapolirane")
plot!(x1,y1,label="Stvarne")

itp = scale(interpolate(ydata2,BSpline(Linear()), OnGrid()),LinRange(-1,1,10))
xx2=[LinRange(-3,3,100);];
etpf = extrapolate(itp,Interpolations.Periodic())
yy2=etpf(xx2)
scatter(xdata2,ydata2, label="Početne tačke")
plot!(xx2,yy2,label="Ekstrapolirane")
plot!(x2,y2,label="Stvarne")

itp = scale(interpolate(ydata1,BSpline(Linear()), OnGrid()),xdata1[1]:(xdata1[2]-xdata1[1]):xdata1[end])
xx1=[LinRange(-1pi,3pi,100);];
etpf = extrapolate(itp,Interpolations.Reflect())
yy1=etpf(xx1)
scatter(xdata1,ydata1, label="Početne tačke")
plot!(xx1,yy1,label="Ekstrapolirane")
plot!(x1,y1,label="Stvarne")

itp = scale(interpolate(ydata2,BSpline(Linear()), OnGrid()),LinRange(-1,1,10))
xx2=[LinRange(-3,3,100);];
etpf = extrapolate(itp,Interpolations.Reflect())
yy2=etpf(xx2)
scatter(xdata2,ydata2, label="Početne tačke")
plot!(xx2,yy2,label="Ekstrapolirane")
plot!(x2,y2,label="Stvarne")

itp = scale(interpolate(ydata1,BSpline(Linear()), OnGrid()),xdata1[1]:(xdata1[2]-xdata1[1]):xdata1[end])
xx1=[LinRange(-1pi,3pi,100);];
etpf = extrapolate(itp,8)
yy1=etpf(xx1)
scatter(xdata1,ydata1, label="Početne tačke")
plot!(xx1,yy1,label="Ekstrapolirane")
plot!(x1,y1,label="Stvarne")

itp = scale(interpolate(ydata2,BSpline(Linear()), OnGrid()),LinRange(-1,1,10))
xx2=[LinRange(-3,3,100);];
etpf = extrapolate(itp,5)
yy2=etpf(xx2)
scatter(xdata2,ydata2, label="Početne tačke")
plot!(xx2,yy2,label="Ekstrapolirane")
plot!(x2,y2,label="Stvarne")