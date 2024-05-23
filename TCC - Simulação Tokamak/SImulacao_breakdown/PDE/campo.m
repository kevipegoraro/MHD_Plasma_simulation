function out=campo(R,Z)
load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/green_table/green_table_nova_65x65_2.mat');

I1 = 1*100; %corrente passando pelas bobinas
fac = 0.4; % fator de mudança da corrente das bobinas verticias internas para as bobinas verticias externas.
% Current moments
Br   = zeros(size(E1.G)); % definindo os Br e bz como zero
Bz = Br;
R0 = 0.37;
B0 = 0.7; 
Vloop = 10;
[rg,~] = meshgrid(r,z);  % criando a mash grid do r e z.
Iq   = [-I1 I1 I1*fac -I1*fac]*1e4; % Quadrupole field
% calculando os campos das 4 bobinas verticais primarias 
for icoil = 1:4
    eval(['Br = Br + E' num2str(icoil) '.BR*Iq(' num2str(icoil) ');'])
    eval(['Bz = Bz + E' num2str(icoil) '.BZ*Iq(' num2str(icoil) ');'])
end
Bphi = R0*B0./rg;
Ephi = -Vloop/2/pi./rg;
%keyboard
out.Br = interp2(r,z,Br,R,Z,'cubic');
out.Bz = interp2(r,z,Bz,R,Z,'cubic');
out.Bphi = interp2(r,z,Bphi,R,Z,'cubic');
out.Ephi = interp2(r,z,Ephi,R,Z,'cubic');