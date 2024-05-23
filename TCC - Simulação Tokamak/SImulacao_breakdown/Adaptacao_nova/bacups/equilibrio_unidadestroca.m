
%build_greentable_nova criar novamente a tabela de campos
%build_greentable_nova
load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/green_table/green_table_nova_65x65.mat')
cnames = {'E1' 'E2' 'E3' 'E4' 'Ep0' 'Ep1' 'Ep2' 'Ep3' 'Ep4' 'E6' 'E7' 'E8'  'E9' 'E10'};

lines = 65;% numero de linhas no contour
ip=-500*1/8  %corrente plasma
iaux_1=500*0.5 %bobina de compensação externa
iaux_2=500*0.8 %bobina de compensação interna
i1=500*2
i2=i1/2 %bobinas verticais
pp = (4*Ep0.G+Ep1.G+Ep2.G+Ep3.G+Ep4.G)*ip; %psi plasma
paux_1 = (E5.G+E6.G+E7.G+E8.G)*iaux_1;
paux_2 = (E9.G+E10.G)*iaux_2;
pv = (E1.G+E2.G)*i1-(E3.G+E4.G)*i2; %psi vacuo


colormap(jet);
figure(1); 
clf; 
contour(r,z,(pp+pv+paux_1+paux_2).',lines); 
axis equal; 
plot_nova(1); 
%kevi = num2str(iaux_2);
%text(0.43,0.45,kevi);
%text(0,0.45,'corrente b.comenp. int.');
colorbar;

%%
%salvar
Campos = (pp+pv+paux_1+paux_2).';
fname = ['/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/green_table/Campos.mat'];
save('Campos.mat')

%%


% bzp = (Ep0.BZ+Ep1.BZ+Ep2.BZ+Ep3.BZ+Ep4.BZ)*ip; %psi plasma
% 
% bzv = (E1.BZ+E2.BZ)*i1-(E3.BZ+E4.BZ)*i2; %psi vacuo
% bzaux_1 = (E6.BZ+E7.BZ)*iaux_1;
% bzaux_2 = (E8.BZ+E9.BZ)*iaux_2;
% colormap(jet);
% figure(1); 
% clf; 
% T = (pp+pv+paux_1+paux_2).';
% T2 = (bzp+bzv+bzaux_1+bzaux_2).'/1e-3;
