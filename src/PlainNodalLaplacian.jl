export getNodalLaplacianMatrix
function dxxMat(n::Int64,h::Float64)
BC = 1.0; # one is for first order Neumann, and 2 is for 2nd order. Note that with 2, the matrix is not symmetric!!!
O1 = -ones(n-1);
O1[n-1] = -BC;
O2 = 2.0*ones(n);
O2[1] = BC;
O2[n] = BC;
O3 = -ones(n-1);
O3[1] = -BC;
dxx = spdiagm((O1/(h^2),O2/(h^2),O3/(h^2)),[-1,0,1],n,n);
return dxx;
end

function getNodalLaplacianMatrix(Msh::RegularMesh)
nodes = Msh.n+1;
I1 = speye(nodes[1]);
Dxx1 = dxxMat(nodes[1],Msh.h[1]);
I2 = speye(nodes[2]);
Dxx2 = dxxMat(nodes[2],Msh.h[2]);
if Msh.dim==2
	L = kron(I2,Dxx1) + kron(Dxx2,I1);
else
	I3 = speye(nodes[3]);
	Dxx3 = dxxMat(nodes[3],Msh.h[3]);
	L = kron(I3,kron(I2,Dxx1) + kron(Dxx2,I1)) + kron(Dxx3,kron(I2,I1));
end
return L;
end	



