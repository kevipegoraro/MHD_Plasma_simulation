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
lines = 500;% numero de linhas no contour
ip=500*4/8;  %corrente plasma
iaux_1=500; %bobina de compensação externa
iaux_2=500; %bobina de compensação interna
j=2;
i1=500*j*0.8; i2=i1*0.5; %bobinas verticais
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
clf;  psi=0*pp+pv+(paux_1+paux_2);
contour(r,z,psi.',lines); 
hold on; 
% aux = 1; 
% for i=1:300
%     for j=1:300
%         if abs(psi(i,j))<aux 
%             aux = psi(i,j);
%         end
%     end
% end
    psi_xp = 0.0003528;%max(max(-abs(psi)));     
contour(r,z,psi.',[1 1]*psi_xp,'k','linewidth',1)
    psi_xp = 0.0003582;%max(max(-abs(psi)));     
%contour(r,z,psi.',[1 1]*psi_xp,'k','linewidth',1)
hold off
axis equal; 
plot_nova(1); 
% kevi = num2str(iaux_2);
% text(0.43,0.45,kevi);
% text(0,0.45,'corrente b.comenp. int.');
colorbar;
%pause(1);
hold off
colormap(jet);
figure(2); 
clf;  psi=pp+pv+(paux_1+paux_2);   psi_xp =0.001302;%max(max(-abs(psi))); 
contour(r,z,psi.',lines); axis equal;  %bzp+bzv+bzaux_1+bzaux_2
%hold on; 
%contour(r,z,psi.',[1 1]*psi_xp,'k','linewidth',2)
%hold off
plot_nova(2); %R0+0.005
colorbar;
%hold off;
 g = log(8*0.37/0.05)+0.05+1.2/2-1.5;
 bz = 1e-4*ip*g/0.37;
%text(0.43,0.45,num2str(bz));
%text(0,0.45,'bz nessesário:');
%axis([-1.1 -0.9]*bz);
%pause(1);
%end
%end
%%