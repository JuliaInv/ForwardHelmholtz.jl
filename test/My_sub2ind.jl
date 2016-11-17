
function My_sub2ind(n::Array{Int64},sub::Array{Int64})
if length(sub)==2
	return sub2ind(tuple(n...),sub[1],sub[2]);
else
	return sub2ind(tuple(n...),sub[1],sub[2],sub[3]);
end
end