function [E_x, E_z] = teste_taylor
%syms x
%g = carga/((P(1)^2+P(2)^2+P(3)^2)+a^2-2*a*(P(1)*cos(x)+P(2)*sin(x)));
%t = taylor(g, 'ExpansionPoint', 0, 'Order', 12);
%size(char(t));
%t = simplify(t);
%size(char(t));
%xd = 0:0.01:2*pi;
%yd = subs(g,x,xd);
%fplot(t, [0, 2*pi])
hold on
%plot(xd, yd, 'r-.')
%title('Taylor approximation vs. actual function')
%legend('Taylor','Function')
%fun = @(x,P(1),P(2),P(3),a,carga) carga/((P(1)^2+P(2)^2+P(3)^2)+a^2-2*a*(P(1)*cos(x)+P(2)*sin(x)));
%f = @(x) t
%posições tokamak
r = 0.06;
pos_bv_in = [0.161 0.265];
pos_bv_ex = [0.53 0.34];

p_central = 0.37;
x_0=p_central-r; %posições ponto P arbtrario
y_0=0;
z_0=0;
a = pos_bv_in(1); 
a2 = pos_bv_ex(1); %raio anel
P = [x_0 y_0 z_0];
H=pos_bv_in(2);
H2 = pos_bv_ex(2);
Ng = 40;
t_contour=1;
out.x = linspace(p_central-t_contour*r,p_central+t_contour*r,Ng);
out.z = linspace(-t_contour*r,t_contour*r,Ng);
[X,Z] = meshgrid(out.x,out.z);
Mo = pi*4E-7;
corrente_bobinas_ex = 0;%Mo*4;
corrente_bobinas_in = 0;%-Mo*4; %proporção 1.369
corrente_plasma = Mo*10;
xi = linspace(p_central-t_contour*r,p_central+t_contour*r,Ng); %x_0:0.01:x_0+2*r;
zi =  linspace(-t_contour*r,t_contour*r,Ng); %-r:0.01:r;
E_x = zeros(length(xi),length(zi));
E_z = E_x;
%length(xi);  size(E_x);  length(zi);  size(E_z);
for ii=1:length(xi)
    for jj=1:length(zi)
        P = [xi(ii),0,zi(jj)];
        aux = campo_anel(P,corrente_bobinas_in,a,H)+campo_anel(P,corrente_bobinas_in,a,-H)+campo_anel(P,corrente_bobinas_ex,a2,H2)+campo_anel(P,corrente_bobinas_ex,a2,-H2);
        E =  corrente_plasma/((xi(ii)-p_central)^2+(zi(jj))^2);
        raio = sqrt((xi(ii)-p_central)^2+(zi(jj))^2);
        campo_corrente_plasma_x = E*((xi(ii)-p_central)/(raio+0.001));
        campo_corrente_plasma_z = E*(zi(jj)/(raio+0.001));
        E_x(ii,jj) = aux(1)+campo_corrente_plasma_x;
        E_z(ii,jj) = aux(3)+campo_corrente_plasma_z;
    end
end
%[DX,DZ] = gradient(B,.24,.24,.24);
quiver(X,Z,E_x,E_z,5)
%contour(X,Z,[E_x,E_z])
end
