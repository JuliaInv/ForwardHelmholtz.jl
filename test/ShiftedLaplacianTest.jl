using jInv.Mesh;
using ForwardHelmholtz
using Multigrid

include("My_sub2ind.jl");

plotting = false;

if plotting
	using PyPlot;
	close("all");
end

# m = readdlm("SEGmodel2Dsalt.dat"); m = m'; m = m*1e-3;
m = 1.5*ones(257,129);


Minv = getRegularMesh([0.0,13.5,0.0,4.2],collect(size(m))-1);

pad = 16;


m = 1./m.^2

f = 2.5;

w = 2*pi*f

println("omega*h:");
println(w*Minv.h*sqrt(maximum(m)));
pad = pad*ones(Int64,Minv.dim);

maxOmega = getMaximalFrequency(m,Minv);
ABLamp = maxOmega;


Sommerfeld = true;
NeumannAtFirstDim = true;
H,gamma = GetHelmholtzOperator(Minv,m,w,w*ones(size(m))*0.01,NeumannAtFirstDim,pad,ABLamp,Sommerfeld);
shift = [0.2;0.2;0.2];
SH = H + GetHelmholtzShiftOP(m, real(w),shift[1]);



n = Minv.n+1; n_tup = tuple(n...);
src = div.(n,2);
src[end] = 1;
q = zeros(Complex,n_tup)
q[My_sub2ind(n,src)] = 1/(Minv.h[1]^2);

levels      = 3;
numCores 	= 2; 
maxIter     = 30;
relativeTol = 1e-6;
# relaxType   = "Jac-GMRES";
relaxType   = "Jac";
relaxParam  = 0.75;
relaxPre 	= 2;
relaxPost   = 2;
cycleType   ='W';
coarseSolveType = "NoMUMPS";

MG = getMGparam(levels,numCores,maxIter,relativeTol,relaxType,relaxParam,
				relaxPre,relaxPost,cycleType,coarseSolveType,0.5,0.0)


				
## Preparing a point shource RHS ############################
n = Minv.n+1; n_tup = tuple(n...);
src = div.(n,2);
src[end] = 1;
q = zeros(Complex128,n_tup)
q[My_sub2ind(n,src)] = 1/(Minv.h[1]^2);
###############################################
b = q[:];

Hparam = HelmholtzParam(Minv,gamma,vec(m),w,NeumannAtFirstDim,Sommerfeld);
Ainv = getShiftedLaplacianMultigridSolver(Hparam, MG,shift,"GMRES",5,true);
Ainv = copySolver(Ainv);

tic()
x = solveLinearSystem(SH',b,Ainv)[1];
toc()
MG.relaxType = "Jac";
MG.cycleType = 'W';
Ainv = getShiftedLaplacianMultigridSolver(Hparam, MG,shift,"BiCGSTAB",0,true);

tic()
x = solveLinearSystem(SH',b,Ainv)[1];
toc()

if plotting
	reX = real(reshape(x,n_tup));
	reX[My_sub2ind(n,src)] = 0.0;
	println(norm(q[:]-H*x)/norm(q[:]));
	println(norm(q[:]-H*y)/norm(q[:]));
	figure();
	imshow(reX');title("Helmholtz Iterative Solution");
end

if plotting
	s = H\q[:];
	s = real(reshape(s,n_tup));
	s[My_sub2ind(n,src)] = 0.0;
	figure();
	imshow(s');title("Helmholtz True Solution");
end

# println("Doing Transpose")
# y = solveLinearSystem(SH',b,Ainv,1)[1];

# println(vecnorm(H'*y-b)/vecnorm(b));


################## WITHOUT THE COARSEST GRID SOL ########################
# coarseSolveType = "DDNodal";
# MG = getMGparam(levels,numCores,maxIter,relativeTol,relaxType,relaxParam,
				# relaxPre,relaxPost,cycleType,coarseSolveType,0.5,0.0)
# Ainv = getShiftedLaplacianMultigridSolver(Hparam, MG,shift,"BiCGSTAB",0,true);


# y = solveLinearSystem(SH',b,Ainv)[1];
###########################################################################


println("SOLVING MULTIPLE RHSs")
coarseSolveType = "Julia";
MG = getMGparam(levels,numCores,maxIter,relativeTol,relaxType,relaxParam,
				relaxPre,relaxPost,cycleType,coarseSolveType,0.5,0.0)
Ainv = getShiftedLaplacianMultigridSolver(Hparam, MG,shift,"BiCGSTAB",0,true);
nrhs = 2;
b = rand(Complex128,length(b),nrhs);
clear!(Ainv);
Ainv.helmParam = Hparam;
tic()
x = solveLinearSystem(SH',b,Ainv)[1];
toc()
println(vecnorm(H*x - b)/vecnorm(b))

MG = getMGparam(levels,numCores,maxIter,relativeTol,relaxType,relaxParam,relaxPre,relaxPost,cycleType,coarseSolveType);
MG.relaxType = "Jac-GMRES";
MG.cycleType = 'K';
Ainv = getShiftedLaplacianMultigridSolver(Hparam, MG,shift,"GMRES",5,true);
tic()
x = solveLinearSystem(SH',b,Ainv)[1];
toc()
println(vecnorm(H*x - b)/vecnorm(b))
