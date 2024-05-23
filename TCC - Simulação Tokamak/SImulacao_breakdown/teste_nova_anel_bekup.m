%testes com campo da bobina sendo um anel
function teste_nova
close all
%%
%load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/green_table/green_table_nova_300x300.mat');

%load('Campos.mat');
%%
Ng = 10;
R1 = 0.06;
Z0 = 0;
dZ = 0.26;
NR = 2;
r0=0.12/2; 
%posições tokamak
r = 0.06;
pos_bv_in = [0.161 0.265];
pos_bv_ex = [0.53 0.34];
R_ini = 0; 
R_fim = 0.6;
Z_ini = -(R_fim-R_ini); 
Z_fim = (R_fim-R_ini); 
p_central = 0.37;

t_contour=2;
out.r = linspace(-p_central-t_contour*r,p_central+t_contour*r,Ng);
out.y = linspace(-p_central-t_contour*r,p_central+t_contour*r,Ng);
out.z = linspace(-t_contour*r,t_contour*r,Ng);
[R,Y,Z] = meshgrid(out.r,out.y,out.z);
raio = linspace(0.001,2*r0,Ng);
theta = linspace(0,2*pi,Ng);
fhi = linspace(0,2*pi,Ng);
%cria a mashgrid psedoridal
%[b11, b12, b13] = meshgrid(raio,theta, fhi);
%[R,Z] = meshgrid(out.r,out.z);
Mo = pi*4E-7;
corrente_bobinas_ex = Mo*4;
corrente_bobinas_in = Mo*2.92; %proporção 1.369
corrente_plasma = Mo*0.07;
%calculando cada um dos N campos da bobina pos_bv_in
N=24;
%B inicial angulo=0
r_aux = 0; %pos_bv_in(1)*sin(0);
y_aux = 0; %pos_bv_in(1)*cos(0);
r_aux2 = 0; %pos_bv_ex(1)*sin(0);
y_aux2 = 0; %pos_bv_ex(1)*cos(0);
%B= corrente_bobinas_in*(1./sqrt((R-r_aux).^2+(Z+pos_bv_in(2)).^2+(Y-y_aux).^2)+1./sqrt((R-r_aux).^2+(Z-pos_bv_in(2)).^2+(Y-y_aux).^2))+corrente_bobinas_ex*(1./sqrt((R-r_aux2).^2+(Z+pos_bv_in(2)).^2+(Y-y_aux2).^2)+1./sqrt((R-r_aux2).^2+(Z-pos_bv_ex(2)).^2+(Y-y_aux2).^2))+corrente_plasma./sqrt((R-p_central*sin(0)).^2+(Z).^2+(Y-p_central*cos(0)).^2);
% for  i=1:1:N %for que ira passar por todos os pedaços dq da bobina.
%     angulo = 2*pi*(i/N);
%     r_aux = pos_bv_in(1)*sin(angulo);
%     y_aux = pos_bv_in(1)*cos(angulo);
%     r_aux2 = pos_bv_ex(1)*sin(angulo);
%     y_aux2 = pos_bv_ex(1)*cos(angulo);
%     F_plasma=corrente_plasma./sqrt((R-p_central*sin(angulo)).^2+(Z).^2+(Y-p_central*cos(angulo)).^2);
%     aux = corrente_bobinas_in*(1./sqrt((R-r_aux).^2+(Z+pos_bv_in(2)).^2+(Y-y_aux).^2)+1./sqrt((R-r_aux).^2+(Z-pos_bv_in(2)).^2+(Y-y_aux).^2))+corrente_bobinas_ex*(1./sqrt((R-r_aux2).^2+(Z+pos_bv_in(2)).^2+(Y-y_aux2).^2)+1./sqrt((R-r_aux2).^2+(Z-pos_bv_ex(2)).^2+(Y-y_aux2).^2));
%     B=B+aux+F_plasma;
% end
    B= corrente_bobinas_in*(1./sqrt((R-r_aux).^2+(Z+pos_bv_in(2)).^2+(Y-y_aux).^2)+1./sqrt((R-r_aux).^2+(Z-pos_bv_in(2)).^2+(Y-y_aux).^2));
   
    [DX,DY,DZ] = gradient(B,.024,.024,.024);
    %keyboard
    quiver3(R,Y,Z,DX,DY,DZ)
   
    
axis equal
grid on
axis([-R_fim R_fim -R_fim R_fim Z_ini Z_fim])
%contour(R,Z,F_plasma,50)
xlabel('R (m)')
ylabel('Z (m)')
hold on
%plots


alfa = 0:0.2:(2*pi+0.1);
xx = cos(alfa)*r;
zz = sin(alfa)*r;
r_poloidal = 0.075;
xx2 = cos(alfa)*r_poloidal;
zz2 = sin(alfa)*r_poloidal;
for i=0:1:N+1
    angulo = 2*pi*(i/N);
    plot3(pos_bv_in(1)*sin(angulo),pos_bv_in(1)*cos(angulo),pos_bv_in(2),'xr','linewidth',3)
    plot3(pos_bv_in(1)*sin(angulo),pos_bv_in(1)*cos(angulo),-pos_bv_in(2),'xr','linewidth',3)
    plot3(pos_bv_ex(1)*sin(angulo),pos_bv_ex(1)*cos(angulo),pos_bv_ex(2),'xr','linewidth',3)
    plot3(pos_bv_ex(1)*sin(angulo),pos_bv_ex(1)*cos(angulo),-pos_bv_ex(2),'xr','linewidth',3)
    r_aux = p_central*sin(angulo);
    y_aux = p_central*cos(angulo);
    plot3(r_aux+xx*sin(angulo),ones(1,length(alfa))*y_aux+xx*cos(angulo),zz,'g','linewidth',2)
    plot3(r_aux+xx2*sin(angulo),ones(1,length(alfa))*y_aux+xx2*cos(angulo),zz2,'b','linewidth',1)
end
%keyboard
%%
 %[DX,DY] = gradient(Campos,.024,.024);
    %keyboard
 %quiver(R,Y,DX,DY) 
 %quiver3(R,Y,Z,DX,DY,DZ)
 %%
end
