%testes com campo da bobina sendo um anel
function teste_nova_anel
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
R_ini = 0; 
R_fim = 0.6;
Z_ini = -(R_fim-R_ini); 
Z_fim = (R_fim-R_ini); 
p_central = 0.37;

t_contour=1;
out.r = linspace(-p_central-t_contour*r,p_central+t_contour*r,Ng);
out.y = linspace(-p_central-t_contour*r,p_central+t_contour*r,Ng);
out.z = linspace(-t_contour*r,t_contour*r,Ng/10);
[R,Y,Z] = meshgrid(out.r,out.y,out.z);
%[R,Z] = meshgrid(out.r,out.z);
Mo = pi*4E-7;
corrente_bobinas_ex = Mo*4;
corrente_bobinas_in = -Mo*2.92; %proporção 1.369
corrente_plasma = Mo*7;
%calculando cada um dos N campos da bobina pos_bv_in
N=48;
%B inicial angulo=0
r_aux = pos_bv_in(1)*sin(0);
y_aux = pos_bv_in(1)*cos(0);
r_aux2 = pos_bv_ex(1)*sin(0);
y_aux2 = pos_bv_ex(1)*cos(0);
B= corrente_bobinas_in*(1./sqrt((R-r_aux).^2+(Z+pos_bv_in(2)).^2+(Y-y_aux).^2)+1./sqrt((R-r_aux).^2+(Z-pos_bv_in(2)).^2+(Y-y_aux).^2))+corrente_bobinas_ex*(1./sqrt((R-r_aux2).^2+(Z+pos_bv_in(2)).^2+(Y-y_aux2).^2)+1./sqrt((R-r_aux2).^2+(Z-pos_bv_ex(2)).^2+(Y-y_aux2).^2))+corrente_plasma./sqrt((R-p_central*sin(0)).^2+(Z).^2+(Y-p_central*cos(0)).^2);
%ZZ = corrente_bobinas_in*(1./sqrt((R-r_aux).^2+(Y-y_aux).^2)+1./sqrt((R-r_aux).^2+(Y-y_aux).^2))+corrente_bobinas_ex*(1./sqrt((R-r_aux2).^2+(Y-y_aux2).^2)+1./sqrt((R-r_aux2).^2+(Y-y_aux2).^2))+corrente_plasma./sqrt((R-p_central*sin(0)).^2+(Y-p_central*cos(0)).^2);
for  i=1:1:N %for que ira passar por todos os pedaços dq da bobina.
    angulo = 2*pi*(i/N);
    r_aux = pos_bv_in(1)*sin(angulo);
    y_aux = pos_bv_in(1)*cos(angulo);
    r_aux2 = pos_bv_ex(1)*sin(angulo);
    y_aux2 = pos_bv_ex(1)*cos(angulo);
    F_plasma=corrente_plasma./sqrt((R-p_central*sin(angulo)).^2+(Z).^2+(Y-p_central*cos(angulo)).^2);
    aux = corrente_bobinas_in*(1./sqrt((R-r_aux).^2+(Z+pos_bv_in(2)).^2+(Y-y_aux).^2)+1./sqrt((R-r_aux).^2+(Z-pos_bv_in(2)).^2+(Y-y_aux).^2))+corrente_bobinas_ex*(1./sqrt((R-r_aux2).^2+(Z+pos_bv_in(2)).^2+(Y-y_aux2).^2)+1./sqrt((R-r_aux2).^2+(Z-pos_bv_ex(2)).^2+(Y-y_aux2).^2));
    B=B+aux+F_plasma;
    %ZZ= ZZ+corrente_bobinas_in*(1./sqrt((R-r_aux).^2+(Y-y_aux).^2)+1./sqrt((R-r_aux).^2+(Y-y_aux).^2))+corrente_bobinas_ex*(1./sqrt((R-r_aux2).^2+(Y-y_aux2).^2)+1./sqrt((R-r_aux2).^2+(Y-y_aux2).^2))+corrente_plasma./sqrt((R-p_central*sin(angulo)).^2+(Y-p_central*cos(angulo)).^2);
end
    [DX,DY,DZ] = gradient(B,.24,.24,.24);
    quiver3(R,Y,Z,DX,DY,DZ)
axis equal
grid on
axis([-R_fim R_fim -R_fim R_fim Z_ini Z_fim])
%contour3(ZZ)
xlabel('R (m)')
ylabel('Z (m)')
hold on
%plots
ploto=1;
if ploto == 1
alfa = 0:0.2:(2*pi+0.1);
xx = cos(alfa)*r;
zz = sin(alfa)*r;
r_poloidal = 0.075;
xx2 = cos(alfa)*r_poloidal;
zz2 = sin(alfa)*r_poloidal;
%N=2*N;
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
end

end
