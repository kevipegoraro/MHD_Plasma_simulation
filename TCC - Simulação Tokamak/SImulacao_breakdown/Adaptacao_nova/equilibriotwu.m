%%
%build_greentable_nova criar novamente a tabela de campos
%build_greentable_nova
% Loading Green's functions table and plasma poloidal flux distribution
load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/green_table/green_table_nova_300x300.mat');
%%
%load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/green_table/green_table_nova_65x65.mat');

%B_t - campo magnético toroidal maximo 15kG
%for j = 3:2:10
%for t = 3:0.5:5
lines = 300;% numero de linhas no contour
ip=500*13.5/8*0;  %corrente plasma
iaux_1=500*3.5; %bobina de compensação externa
iaux_2=500*4; %bobina de compensação interna
j=2;
i1=500*5*j; i2=i1*1.1; %bobinas verticais
pp = (4*E50.G+E51.G+E52.G+E53.G+E54.G)*ip; %psi plasma
paux_1 = (E6.G+E7.G)*iaux_1;
paux_2 = (E8.G+E9.G)*iaux_2;
pv = (E1.G+E2.G)*i1-(E3.G+E4.G)*i2; %psi vacuo

bzp = (4*E50.BZ+E51.BZ+E52.BZ+E53.BZ+E54.BZ)*ip; %psi plasma

bzv = (E1.BZ+E2.BZ)*i1-(E3.BZ+E4.BZ)*i2; %psi vacuo
bzaux_1 = (E6.BZ+E7.BZ)*iaux_1;
bzaux_2 = (E8.BZ+E9.BZ)*iaux_2;
%%
colormap(jet);
figure(1); 
clf; 
contour(r,z,(pp+pv+paux_1+paux_2).',lines); 
axis equal; 
plot_nova(1); 

colormap(jet);
figure(3); 
clf; 
ka=(4*E50.BR+E51.BR+E52.BR+E53.BR+E54.BR)*ip;
kaa=(4*E50.BZ+E51.BZ+E52.BZ+E53.BZ+E54.BZ)*ip;
contour(r,z,sqrt(ka.^2+kaa.^2).'); 
axis equal; 
plot_nova(3); 
% kevi = num2str(iaux_2);
% text(0.43,0.45,kevi);
% text(0,0.45,'corrente b.comenp. int.');
colorbar;
%pause(1);
hold off
colormap(jet);
figure(2); 
clf; 
contourf(r,z,(bzp+bzv+bzaux_1+bzaux_2).'/1e-3,lines); axis equal; 
plot_nova(2); 
colorbar;
%hold off;
 g = log(8*0.37/0.05)+0.05+1.2/2-1.5;
 bz = 1e-4*ip*g/0.37;
text(0.43,0.45,num2str(bz));
text(0,0.45,'bz nessesário:');
caxis([-1.1 -0.9]*bz);
%pause(1);
%end
%end
%%