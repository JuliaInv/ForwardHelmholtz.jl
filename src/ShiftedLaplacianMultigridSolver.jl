export ShiftedLaplacianMultigridSolver,updateParam,getShiftedLaplacianMultigridSolver,copySolver
export solveLinearSystem,clear!

type ShiftedLaplacianMultigridSolver <: AbstractSolver
	helmParam		:: HelmholtzParam
	MG				:: MGparam
	shift			:: Array{Float64}
	Krylov			:: String
	inner			:: Int64
	doClear			:: Int64 # flag to clear factorization
	verbose			:: Bool
	setupTime::Real
	nPrec::Int
	solveTime::Real
end

import jInv.LinearSolvers.copySolver;
function copySolver(s::ShiftedLaplacianMultigridSolver)
	# copies absolutely what's necessary.
	clear!(s.helmParam.Mesh);
	return getShiftedLaplacianMultigridSolver(s.helmParam,Multigrid.copySolver(s.MG),s.shift,s.Krylov,s.inner,s.verbose);
end

function getShiftedLaplacianMultigridSolver(helmParam::HelmholtzParam, MG::MGparam,shift::Array{Float64},Krylov::String="BiCGSTAB",inner::Int64=5,verbose::Bool = false)
	return ShiftedLaplacianMultigridSolver(helmParam,MG,shift,Krylov,inner,0,verbose,0.0,0,0.0);
end

function getShiftedLaplacianMultigridSolver(helmParam::HelmholtzParam, MG::MGparam,shift::Float64,Krylov::String="BiCGSTAB",inner::Int64=5,verbose::Bool = false)
	return getShiftedLaplacianMultigridSolver(helmParam,MG,ones(MG.levels)*shift,Krylov,inner,verbose);
end

function solveLinearSystem(ShiftedHT,B,param::ShiftedLaplacianMultigridSolver,doTranspose::Int64=0)
	if size(B,2) == 1
		B = vec(B);
	end
	if param.doClear==1
		clear!(param.MG);
	end
	if vecnorm(B) == 0.0
		X = zeros(eltype(B),size(B));
		return X, param;
	end
	tic()
	TYPE = eltype(B);
	n = size(B,1)
	nrhs = size(B,2)

	# build preconditioner
	if hierarchyExists(param.MG)==false
		ShiftedLaplacianMGsetup(GetHelmholtzOperator, param.helmParam,param.MG,param.shift,TYPE,nrhs,param.verbose);
		# MGsetup(ShiftedHT,param.helmParam.Mesh,param.MG,TYPE,nrhs,param.verbose);
	end

	if (doTranspose != param.MG.doTranspose)
		transposeHierarchy(param.MG);
	end
	
	adjustMemoryForNumRHS(param.MG,TYPE,size(B,2),param.verbose)
	BLAS.set_num_threads(param.MG.numCores);
	ShiftedHT = param.MG.As[1];
	Az = param.MG.memCycle[1].b;
	if param.shift[1] != 0.0
		MShift = GetHelmholtzShiftOP(param.helmParam.m, param.helmParam.omega,param.shift[1]);
		# MShift = GetHelmholtzShiftOPNew(param.helmParam.m, param.helmParam.omega,param.shift[1]);
		if doTranspose==1
			MShift = -MShift; # this is because of the conjugate
			Hfun = getHelmholtzFun(ShiftedHT,MShift,param.MG.numCores);
			# println("doTranspose")
		else
			Hfun = getHelmholtzFun(ShiftedHT,MShift,param.MG.numCores);
		end
		function Afun(z::ArrayTypes)
			Hfun(one(eltype(z)),z,zero(eltype(z)),Az);
			return Az;
		end
	else
		Afun = getAfun(ShiftedHT,Az,param.MG.numCores);
	end
	
	param.setupTime += toq();
	tic()
	if param.Krylov=="GMRES"
		X, param.MG,num_iter,nprec = solveGMRES_MG(Afun,param.MG,B,Array(eltype(B),0),param.verbose,param.inner);
	elseif param.Krylov=="BiCGSTAB"
		X, param.MG,num_iter,nprec = solveBiCGSTAB_MG(Afun,param.MG,B,Array(eltype(B),0),param.verbose);
	end
	param.solveTime +=toq();
	param.nPrec += nprec;
	if num_iter >= param.MG.maxOuterIter
		warn("MG solver reached maximum iterations without convergence");
	end

	return X, param
end



import jInv.Utils.clear!
function clear!(s::ShiftedLaplacianMultigridSolver)
	 clear!(s.MG);
	 clear!(s.helmParam)
	 s.doClear = 0;
end

function ShiftedLaplacianMGsetup(A::Function, helmParam::HelmholtzParam,param::MGparam,shift::Array{Float64},rhsType::DataType = Float64,nrhs::Int64 = 1,verbose::Bool=false)
Ps      	= Array(SparseMatrixCSC,param.levels-1);
Rs      	= Array(SparseMatrixCSC,param.levels-1);
As 			= Array(SparseMatrixCSC,param.levels);
Meshes  	= Array(RegularMesh,param.levels); 
relaxPrecs 	= Array(SparseMatrixCSC,param.levels);
n = helmParam.Mesh.n + 1; # n here is the number of NODES!!!
N = prod(n);
As[1] = A(helmParam.Mesh,helmParam.m,helmParam.omega,helmParam.gamma+shift[1]*real(helmParam.omega),helmParam.NeumannOnTop,helmParam.Sommerfeld)';
# As[1] = A(helmParam.Mesh,helmParam.m,helmParam.omega + 1im*0.5*shift[1]*real(helmParam.omega),helmParam.gamma,helmParam.NeumannOnTop,helmParam.Sommerfeld)';
Meshes[1] = helmParam.Mesh;

mNodal = helmParam.m[:];
gamma = helmParam.gamma[:];

Cop = nnz(As[1]);
for l = 1:(param.levels-1)
	if verbose
		tic()
	end
    AT = As[l];
    if param.relaxType=="Jac" || param.relaxType=="Jac-GMRES"
		d = param.relaxParam./diag(AT);
		relaxPrecs[l] = spdiagm(d);# here we need to take the conjugate for the SpMatVec, but we give At instead of A so it cancels
	elseif param.relaxType=="SPAI"
		relaxPrecs[l] = spdiagm(param.relaxParam*getSPAIprec(AT)); # here we need to take the conjugate for the SpMatVec, but we give At instead of A so it cancels
	else
		error("Unknown relaxation type !!!!");
	end
	
	(P,nc) = getFWInterp(n,true);
	if (size(P,1)==size(P,2))
		if verbose; println(string("Stopped Coarsening at level ",l)); end
		param.levels = l;
		As = As[1:l];
		Ps = Ps[1:l-1];
		Rs = Rs[1:l-1];
		Meshes = Meshes[1:l];
		relaxPrecs = relaxPrecs[1:l];
		break;
	else
		RT = (0.5^Meshes[l].dim)*P;
		Rs[l] = RT;# this is becasue we hold the transpose of the matrices and P = R' anyway here....
		Ps[l] = P';
		
		# mNodal = RT'*mNodal[:];
		# gamma = RT'*gamma[:];
		Meshes[l+1] = getRegularMesh(Meshes[l].domain,nc-1);
		# Act = A(Meshes[l+1],mNodal,helmParam.omega,gamma+shift[l+1]*real(helmParam.omega),helmParam.NeumannOnTop,helmParam.Sommerfeld)';
		Act = Ps[l]*AT*Rs[l]; #- GetHelmholtzShiftOP(mNodal, real(helmParam.omega),(shift[l+1] - shift[l])*real(helmParam.omega));
		
		# LAP = getNodalLaplacianMatrix(Meshes[l]);
		# LAPc = getNodalLaplacianMatrix(Meshes[l+1]);
		# Act = Ps[l]*LAP*Rs[l] +  A(Meshes[l+1],mNodal,helmParam.omega,gamma+shift[l+1]*real(helmParam.omega),helmParam.NeumannOnTop,helmParam.Sommerfeld)' - LAPc;
		
		As[l+1] = Act;
		
		Cop = Cop + nnz(Act);
		
		if verbose; println("MG setup: ",n," took:",toq()); end;
		n = nc;
		N = prod(nc);
	end
end
if verbose 
	tic()
	println("MG Setup: Operator complexity = ",Cop/nnz(As[1]));
end
defineCoarsestAinv(param,As[end]);
if verbose 
	println("MG setup coarsest ",param.coarseSolveType,": ",n,", done LU in ",toq());
end
param.As = As;
param.Ps = Ps;
param.Rs = Rs;
param.relaxPrecs = relaxPrecs;
param = adjustMemoryForNumRHS(param,rhsType,nrhs,verbose);
return;
end









