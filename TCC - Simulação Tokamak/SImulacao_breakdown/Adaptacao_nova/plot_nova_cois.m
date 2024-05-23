function plot_nova_cois
%posições tokamak
r = 0.06;
pos_bv_in = [0.161 0.265];
pos_bv_ex = [0.53 0.34];
R_ini = 0; R_fim = 0.53;
Z_ini = -(R_fim-R_ini)/2; %-0.09; 
Z_fim = (R_fim-R_ini)/2; %0.34;
p_central = 0.37;

axis([R_ini R_fim R_ini R_fim Z_ini Z_fim])
%4 bobinas verticais anel
plot(pos_bv_in(1),pos_bv_in(2),'xr','linewidth',3)
plot(pos_bv_in(1),-pos_bv_in(2),'xr','linewidth',3)
plot(pos_bv_ex(1),pos_bv_ex(2),'xr','linewidth',3)
plot(pos_bv_ex(1),-pos_bv_ex(2),'xr','linewidth',3)
%circulo que representa a vazilha de vacuo
alfa = 0:0.1:(2*pi+0.1);
xx = cos(alfa)*r;
yy = sin(alfa)*r;
plot(p_central+xx,yy,'r','linewidth',2)
%bobinas paralelas
r_poloidal = 0.075;
xx = cos(alfa)*r_poloidal;
yy = sin(alfa)*r_poloidal;
plot(p_central+xx,yy,'b','linewidth',2)
plot([p_central+r*(1+1/3) p_central+r*(1+1/3)],[-2*r+0.1 2*r-0.1],'b','linewidth',2)
plot([p_central-r*(3/2) p_central-r*(3/2)],[-2*r+0.1 2*r-0.1],'b','linewidth',2)
end