using LinearAlgebra

A=[2 1 3; 2 6 8; 6 8 18]
B=[1;3;5]
X=A\B
X=(A'*A)\(A'*B)
X=inv(A)*B
X=pinv(A)*B


A=[2 1; -3 1; -1 1]
B=[-1;-1;2]
X=A\B
err=abs.(A*X-B)

A=[2 4 6 8; 3 5 7 9]
B=10*ones(2,4)
A/B
A=[1 2 3; 4 5 6; 7 8 9]
B=[7 5 6; 2 0 8; 5 7 1]
A/B

using Pkg; 
Pkg.add(url="https://github.com/JuliaMatrices/SpecialMatrices.jl")

using SpecialMatrices

C=Riemann(3)
D=Hilbert(3)
det(C)
det(D)
rank(C)
rank(D)
cond(C)
cond(D)


A=[1 -2 -2 -3; 3 -9 0 -9; -1 2 4 7; -3 -6 26 2]
B=[ 2 -1 0; -1 2 -1; 0 -1 2]

L,U,p=lu(A)
A[p,:]
(L*U)
factorize(A)
R=cholesky(B)
R.L*R.U
factorize(B)

M=[3 4 18 34 0 2 31;
1 -3 -7 -6 2 4 26;
2 1 7 16 3 -1 27;
5 11 43 74 2 0 56;
3 -3 -3 6 -1 14 55;
-2 0 -4 -12 1 5 6;
1 -6 -16 -18 4 4 33]

Pkg.add("RowEchelon")
using RowEchelon
rref(M)

P=[1 1 1;1 2 3;1 3 6]
Q,R=qr(P)
Q*R

A=[1 0 1;-1 -2 0;0 1 -1]
X=eigvals(A)
AA=A;
Q,R= qr(AA); AA=R*Q;
Q,R= qr(AA); AA=R*Q;
Q,R= qr(AA); AA=R*Q;
Q,R= qr(AA); AA=R*Q;
Q,R= qr(AA); AA=R*Q;
Q,R= qr(AA); AA=R*Q;
Q,R= qr(AA); AA=R*Q;
Q,R= qr(AA); AA=R*Q;
Q,R= qr(AA); AA=R*Q;
Q,R= qr(AA); AA=R*Q;
Q,R= qr(AA); AA=R*Q;
Q,R= qr(AA); AA=R*Q;
X
AA

F=svd(A)
F.U * Diagonal(F.S) * F.Vt


using Pkg;
Pkg.add("Images")
img_path="lena.png"
Pkg.add("TestImages")
Pkg.add("ImageMagick")

using Images,TestImages,ImageMagick

img = load(img_path)
img1 = Gray.(img)
A = convert(Array{Float64}, img1)
n,m=size(A)
F=svd(A)
U=F.U
S=Diagonal(F.S)
V= F.Vt'
k=255
Ak = U[:,1:k]*S[1:k,1:k]*V[:,1:k]'
Gray.(Ak)