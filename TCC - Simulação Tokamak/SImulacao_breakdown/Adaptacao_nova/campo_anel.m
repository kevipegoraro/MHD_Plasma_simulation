function E =campo_anel(P,carga,a,H)
r=0.06; %6 centimetros de raio
% H  é a altura da bobina vertical
%carga é a carag da bobina vertical
%a é o raio da bobina vertical
%P = [x_0 y_0 z_0]; %posição para calcular campo
%x_0, z_0 vao de -r a r 
P(3)=H-P(3); %correção devido a altura da bobina vertical
ri=(P(1)^2+P(2)^2+P(3)^2)+a^2;
%raio = sqrt(ri-2*a*(P(1)*cos(x)+P(2)*sin(x)));
g_x = @(x)(P(3)/sqrt(ri-2*a*(P(1)*cos(x)+P(2)*sin(x))))*(carga/(ri-2*a*(P(1)*cos(x)+P(2)*sin(x))));
%g_y = @(x)sin(pi/2-(P(3)/(abs(P(2)-a*sin(x))+0.001)))*(carga/(ri-2*a*(P(1)*cos(x)+P(2)*sin(x))));
g_z = @(x)((P(1)-a*cos(x))/sqrt(ri-2*a*(P(1)*cos(x)+P(2)*sin(x))))*(carga/(ri-2*a*(P(1)*cos(x)+P(2)*sin(x))));
x0=0; 
h=0.01; 
l=2*pi;
x = x0:h:l; 
i = x0;
k1=0; 
k2=0; 
k3=0;
for i=1:length(x)
        k1 = k1+g_x(x(i))*h;
        %k2 = k2+g_y(x(i))*h;
        k3 = k3+g_z(x(i))*h;
end
E = [k1 k2 k3];
end

