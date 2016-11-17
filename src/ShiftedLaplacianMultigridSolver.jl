export ShiftedLaplacianMultigridSolver,updateParam,getShiftedLaplacianMultigridSolver,copySolver
export solveLinearSystem,clear!

type ShiftedLaplacianMultigridSolver <: AbstractSolver
	M  				:: RegularMesh
	m				:: Array{Float64}
	omega			:: Float64
	MG				:: MGparam
	shift			:: Float64
	Krylov			:: String
	doClear			:: Int64 # flag to clear factorization
	verbose			:: Bool
end

import jInv.LinearSolvers.copySolver;
function copySolver(s::ShiftedLaplacianMultigridSolver)
	# copies absolutely what's necessary.
	clear!(s.M);
	return getShiftedLaplacianMultigridSolver(s.M,Multigrid.copySolver(s.MG),s.shift,s.Krylov,s.verbose);
end

function getShiftedLaplacianMultigridSolver(M::RegularMesh, MG::MGparam,shift::Float64,Krylov::String="BiCGSTAB",verbose::Bool = false)
	return ShiftedLaplacianMultigridSolver(M,zeros(0),0.0,MG,shift,Krylov,0,verbose);
end

function updateParam(solver::ShiftedLaplacianMultigridSolver,M::RegularMesh,m::Array{Float64},omega:: Float64)
	solver.M     = M;
	solver.m     = m;
	solver.omega = omega;
	return solver;
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
	TYPE = eltype(B);
	n = size(B,1)
	nrhs = size(B,2)

	# build preconditioner
	if hierarchyExists(param.MG)==false
		MGsetup(ShiftedHT,param.MG,TYPE,nrhs,param.verbose);
	end

	if (doTranspose != param.MG.doTranspose)
		transposeHierarchy(param.MG);
	end
	ShiftedHT = param.MG.As[1];
	MShift = GetHelmholtzShiftOP(param.m, param.omega,param.shift);
	BLAS.set_num_threads(param.MG.numCores);
	if doTranspose==1
		MShift = -MShift; # this is because of the conjugate
		Hfun = getHelmholtzFun(ShiftedHT,MShift,param.MG.numCores);
	else
		Hfun = getHelmholtzFun(ShiftedHT,MShift,param.MG.numCores);
	end

	if param.Krylov=="GMRES"
		X, param.MG,num_iter = solveGMRES_MG(Hfun,param.MG,B,Array(eltype(B),0),param.verbose);
	elseif param.Krylov=="BiCGSTAB"
		Az = zeros(eltype(B),size(B));
		function Afun(z::ArrayTypes)
			Hfun(one(eltype(z)),z,zero(eltype(z)),Az);
			return Az;
		end
		X, param.MG,num_iter = solveBiCGSTAB_MG(Afun,param.MG,B,Array(eltype(B),0),param.verbose);
	end

	if num_iter >= param.MG.maxOuterIter
		warn("MG solver reached maximum iterations without convergence");
	end

	return X, param
end
import jInv.Utils.clear!
function clear!(s::ShiftedLaplacianMultigridSolver)
	 Multigrid.clear!(s.MG);
	 clear!(s.M);
	 s.m = zeros(0);
	 s.doClear = 0;
end
