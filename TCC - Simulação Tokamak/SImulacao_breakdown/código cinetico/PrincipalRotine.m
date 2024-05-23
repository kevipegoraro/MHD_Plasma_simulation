function PrincipalRotine()
%%
clear all; clc;
syms f1(v) 
r=1; %[1 2 3 4];
times = 1;
mesh = 1; %[0.37 0.37+0.01 0.37+0.02 0.37+0.03];
%G=DefinProblem(r,mesh, times); % f*(partial/partial_t+v*\nabla+*\nabla_v*(F_ext+Forçalorentz)/(m_\alfa) )
const=1; const2=-1; T=1;

vVal = -10:1:10;
TVal = 1;
f1(v) = const*exp(const2*v^2)/(T^(3/2));
dfv=diff(f1,v); d=dfv(3)
keyboard
sol=solve(dfv(vVal)*f1(vVal)-1, f1(vVal))

%%
%G.mag=DefinGreenTable(G);
%G.pol = GeraPol(G);
G.f=AtualizaDistros(G);
%PlotFunctions(G);
keyboard
end