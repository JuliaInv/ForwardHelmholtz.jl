using jInv.Mesh;
using ForwardHelmholtz
using Multigrid
import EikonalInv.expandModelNearest
import EikonalInv.addAbsorbingLayer

function My_sub2ind(n::Array{Int64},sub::Array{Int64})
if length(sub)==2
	return sub2ind(tuple(n...),sub[1],sub[2]);
else
	return sub2ind(tuple(n...),sub[1],sub[2],sub[3]);
end
end

plotting = false;

if plotting
	using PyPlot;
	close("all");
end

# m = readdlm("SEGmodel2Dsalt.dat"); m = m'; m = m*1e-3;
m = 1.5*ones(256,128);

m = EikonalInv.expandModelNearest(m,[256,128],[128,64]);
Minv = getRegularMesh([0.0,13.5,0.0,4.2],collect(size(m))-1);

pad = 16;
(m,Minv) = EikonalInv.addAbsorbingLayer(m,Minv,pad);

m = 1./m.^2

f = 1.25;
w = 2*pi*f

println("omega*h:");
println(w*Minv.h*sqrt(maximum(m)));
pad = pad*ones(Int64,Minv.dim);
H = GetHelmholtzOperator(Minv,m,w,ones(size(m))*0.0001,true,pad,1.0,true)[1];
shift = 0.15;
SH = H + GetHelmholtzShiftOP(m, w,shift);
n = Minv.n+1; n_tup = tuple(n...);
src = div(n,2);
src[end] = 1;
q = zeros(Complex,n_tup)
q[My_sub2ind(n,src)] = 1/(Minv.h[1]^2);

levels      = 3;
numCores 	= 8; 
maxIter     = 50;
relativeTol = 1e-5;
relaxType   = "SPAI";
relaxParam  = 1.0;
relaxPre 	= 2;
relaxPost   = 2;
cycleType   ='W';
coarseSolveType = "NoMUMPS";

MG = getMGparam(levels,numCores,maxIter,relativeTol,relaxType,relaxParam,
				relaxPre,relaxPost,cycleType,coarseSolveType,0.5,0.0,Minv)

# Ainv = getShiftedLaplacianMultigridSolver(Minv, MG,shift,"GMRES");
Ainv = getShiftedLaplacianMultigridSolver(Minv, MG,shift,"BiCGSTAB");
Ainv = updateParam(Ainv,Minv,vec(m),w);

## Preparing a point shource RHS ############################
n = Minv.n+1; n_tup = tuple(n...);
src = div(n,2);
src[end] = 1;
q = zeros(Complex128,n_tup)
q[My_sub2ind(n,src)] = 1/(Minv.h[1]^2);
###############################################
b = q[:];

tic()
x = solveLinearSystem(SH',b,Ainv)[1];
toc()

reX = real(reshape(x,n_tup));
reX[My_sub2ind(n,src)] = 0.0;

println(norm(q[:]-H*x)/norm(q[:]));

if plotting
	figure();
	imshow(reX');title("Helmholtz Iterative Solution");
end

s = H\q[:];

s = real(reshape(s,n_tup));
s[My_sub2ind(n,src)] = 0.0;

if plotting
	figure();
	imshow(s');title("Helmholtz True Solution");
end

println("Doing Transpose")
y = zeros(Complex128,size(b));
y = solveLinearSystem(speye(3),b,Ainv,1)[1];

println(vecnorm(H'*y-b)/vecnorm(b));


# println("SOLVING MULTIPLE RHSs")
# nrhs = 16;
# b = rand(Complex128,length(b),nrhs);
# HelmholtzForward.clear!(Ainv);
# Ainv = updateParam(Ainv,Minv,vec(m),w);
# tic()
# x = solveLinearSystem(SH',b,Ainv)[1];
# toc()
# println(vecnorm(H*x - b)/vecnorm(b))

# MG = getMGparam(levels,numCores,maxIter,innerIter,relativeTol,relaxType,relaxParam,relaxPre,relaxPost,cycleType,coarseSolveType);
# Ainv = getShiftedLaplacianMultigridSolver(Minv, MG,shift,"GMRES");
# Ainv = updateParam(Ainv,Minv,vec(m),w);
# tic()
# x = solveLinearSystem(SH',b,Ainv)[1];
# toc()
# println(vecnorm(H*x - b)/vecnorm(b))







