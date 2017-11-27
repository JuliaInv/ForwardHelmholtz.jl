
module ForwardHelmholtz

using jInv.Mesh;
using Multigrid;

using jInv.LinearSolvers
import jInv.Utils.clear!
import jInv.LinearSolvers.AbstractSolver
import jInv.LinearSolvers.solveLinearSystem


export HelmholtzParam,getShiftedHelmholtzParam
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

function getShiftedHelmholtzParam(p::HelmholtzParam,s::Float64)
	return HelmholtzParam(p.Mesh,p.gamma + s*real(p.omega), p.m,p.omega,p.NeumannOnTop,p.Sommerfeld);
end


include("GetHelmholtz.jl");

include("PlainNodalLaplacian.jl");
include("ShiftedLaplacianMultigridSolver.jl");
#include("PointSourceADR/PointSourceADR.jl");
#include("Elastic/ElasticHelmholtz.jl");
end