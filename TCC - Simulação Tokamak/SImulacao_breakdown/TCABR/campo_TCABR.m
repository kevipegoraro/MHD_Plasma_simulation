function out=campo_TCABR(R,Z)
global CAMPO G
%if length(R)>1
    
I1 = 1*100; %corrente passando pelas bobinas
fac = 0.35; % fator de mudança da corrente das bobinas verticias internas para as bobinas verticias externas.
% Current moments
Br   = zeros(size(CAMPO.E1.G)); % definindo os Br e bz como zero
Bz = Br;
R0 = 0.37;
B0 = 0.7; 
Vloop = 10;
[rg,~] = meshgrid(CAMPO.r,CAMPO.z);  % criando a mash grid do r e z.
Iq   = [-I1 I1 I1*fac -I1*fac -I1 I1 I1*fac -I1*fac I1*fac]*1e4; % Quadrupole field
% calculando os campos das 4 bobinas verticais primarias 
%keyboard
for icoil = 1:7
    eval(['Br = Br + CAMPO.E' num2str(icoil) '.BR*Iq(' num2str(icoil) ');'])
    eval(['Bz = Bz + CAMPO.E' num2str(icoil) '.BZ*Iq(' num2str(icoil) ');'])
end
%Bplx = CAMPO.Ep0.BR;
%Bply = CAMPO.Ep0.BZ;
Bphi = R0*B0./rg;
Ephi = -Vloop/2/pi./rg;
prop=1;
corrige=0;
%keyboard
%out.Bplx = interp2(CAMPO.r,CAMPO.z,Bplx,R*prop+corrige,Z*prop,'cubic');
%out.Bply = interp2(CAMPO.r,CAMPO.z,Bply,R*prop+corrige,Z*prop,'cubic');
out.Br = interp2(CAMPO.r,CAMPO.z,Br,R*prop+corrige,Z*prop,'cubic');
out.Bz = interp2(CAMPO.r,CAMPO.z,Bz,R*prop+corrige,Z*prop,'cubic');
out.Bphi = interp2(CAMPO.r,CAMPO.z,Bphi,R*prop+corrige,Z*prop,'cubic');
out.Ephi = interp2(CAMPO.r,CAMPO.z,Ephi,R*prop+corrige,Z*prop,'cubic');
