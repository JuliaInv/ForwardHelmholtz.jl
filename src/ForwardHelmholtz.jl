
module ForwardHelmholtz

using jInv.Mesh;
using Multigrid;
using jInv.LinearSolvers
import jInv.Utils.clear!
import jInv.LinearSolvers.AbstractSolver
import jInv.LinearSolvers.solveLinearSystem


export HelmholtzParam
type HelmholtzParam
	Mesh   			:: RegularMesh;
	gamma  			:: Array{Float64};
	m               :: Array{Float64};
	omega			:: Union{Float64,Complex128}
	NeumannOnTop	:: Bool
	Sommerfeld 		:: Bool
end


import jInv.Utils.clear!

function clear!(HP::HelmholtzParam)
	clear!(HP.Mesh);
	m = zeros(0);
	gamma = zeros(0);
	omega = 0.0;
end

include("GetHelmholtz.jl");

include("PlainNodalLaplacian.jl");
include("ShiftedLaplacianMultigridSolver.jl");
end