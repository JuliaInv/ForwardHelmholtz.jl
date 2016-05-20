
module ForwardHelmholtz

using jInv.Mesh;
using Multigrid;
using jInv.LinearSolvers
import jInv.Utils.clear
import jInv.LinearSolvers.AbstractSolver
import jInv.LinearSolvers.solveLinearSystem



type HelmholtzParam
	Mesh   			:: RegularMesh;
	gamma  			:: Array{Float64};
	NeumannOnTop	:: Bool
	Sommerfeld 		:: Bool
end

include("GetHelmholtz.jl");
include("PlainNodalLaplacian.jl");
include("ShiftedLaplacianMultigridSolver.jl");
end