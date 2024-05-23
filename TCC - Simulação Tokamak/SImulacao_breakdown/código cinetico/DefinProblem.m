function G=DefinProblem(r,mesh, times)
G.r=r;
tam=size(mesh);
G.times = times; const=1; const2=-1; %gradient(
%v = syms('v'); T = syms('T'); 
syms f(v,T) 
f(v,T) = const*exp(const2*v)/(T^(3/2));
G.mesh = sym('A%d%d',tam,'positive')
for i=G.r
eval(['G.f' num2str(i) '=G.mesh(i)*G.v*exp(G.times);'])
end
end