%testes com malha
function teste_nova
close all
Ng = 100;
R1 = 0.06;
Z0 = 0;
dZ = 0.26;
NR = 2;
%posições tokamak
r = 0.06;
pos_bv_in = [0.161 0.265];
pos_bv_ex = [0.53 0.34];
R_ini = 0; R_fim = 0.53;
Z_ini = -(R_fim-R_ini)/2; %-0.09; 
Z_fim = (R_fim-R_ini)/2; %0.34;
p_central = 0.37;

t_contour=2;
out.r = linspace(p_central-t_contour*r,p_central+t_contour*r,Ng)
out.z = linspace(-t_contour*r,t_contour*r,Ng)
[R,Z] = meshgrid(out.r,out.z);
clf
%contour(out.r,out.z,Gtot,51,'b')
% = green_em_nova(Ro,Zo,R,Z,mode)
Mo = pi*4E-7;
corrente_bobinas_ex = Mo*4;
corrente_bobinas_in = Mo*2.92; %proporção 1.369
corrente_plasma = Mo*0.07;
B1= corrente_bobinas_ex./sqrt((R-pos_bv_ex(1)).^2+(Z+pos_bv_ex(2)).^2);
B2= corrente_bobinas_ex./sqrt((R-pos_bv_ex(1)).^2+(Z-pos_bv_ex(2)).^2);
B3= corrente_bobinas_in./sqrt((R-pos_bv_in(1)).^2+(Z+pos_bv_in(2)).^2);
B4= corrente_bobinas_in./sqrt((R-pos_bv_in(1)).^2+(Z-pos_bv_in(2)).^2);
%campo eletromagnetico gerado pelo plasma
F_plasma=corrente_plasma./sqrt((R-p_central).^2+(Z).^2);
G =-B1-B2-B3-B4+F_plasma;
%hold off
axis equal
grid on
axis([R_ini R_fim Z_ini Z_fim])
contour(R,Z,G,100)
xlabel('R (m)')
ylabel('Z (m)')
%drawnow
%pause(2)
hold on
plot_nova_cois
end
