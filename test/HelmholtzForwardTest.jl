using jInv.Mesh;
using ForwardHelmholtz
#import EikonalInv.expandModelNearest
#import EikonalInv.addAbsorbingLayer
plotting = false;
if plotting
	using PyPlot;
end


function My_sub2ind(n::Array{Int64},sub::Array{Int64})
if length(sub)==2
	return sub2ind(tuple(n...),sub[1],sub[2]);
else
	return sub2ind(tuple(n...),sub[1],sub[2],sub[3]);
end
end


# m = readdlm("SEGmodel2Dsalt.dat"); m = m';
m = ones(256,128);
m = m*1e-3;
#m = EikonalInv.expandModelNearest(m,[256,128],[128,64]);
Minv = getRegularMesh([0.0,13.5,0.0,4.2],collect(size(m))-1);

pad = 16;
#(m,Minv) = EikonalInv.addAbsorbingLayer(m,Minv,pad);

m = ones(size(m));
m = 1./m.^2

if plotting
	figure();
	imshow(1./sqrt(m'))
end
f = 1.0;
w = 2*pi*f
println("omega*h:");
println(w*Minv.h*sqrt(maximum(m)));
pad = pad*ones(Int64,Minv.dim);
H = GetHelmholtzOperator(Minv,m,w,ones(size(m))*0.01,true,pad,1.0,true)[1];
SH = H + GetHelmholtzShiftOP(m, w,0.1);
n = Minv.n+1; n_tup = tuple(n...);
src = div(n,2);
src[end] = 1;
q = zeros(Complex,n_tup)
q[My_sub2ind(n,src)] = 1/(Minv.h[1]^2);

s = H\q[:];
s = real(reshape(s,n_tup));
s[My_sub2ind(n,src)] = 0.0;

if plotting
	figure();
	imshow(s');title("Helmholtz Solution");
end


sh = SH\q[:];
sh = real(reshape(sh,n_tup));
sh[My_sub2ind(n,src)]= 0.0;

if plotting
	figure();
	imshow(sh');title("Shifted Helmholtz Solution");
end