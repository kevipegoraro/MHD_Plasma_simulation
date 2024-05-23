function y1 = RK4(x0,y0,h,l)
x = x0:h:l; 
y1 = zeros(4, length(x));
%y_1 = zeros(length(x));
y1(4,1) = y0(4);
%y_2 = zeros(length(x));
y1(3,1) = y0(3);
%y_3 = zeros(length(x));
y1(2,1) = y0(2);
%y_4 = zeros(length(x));
y1(1,1) = y0(1);
%y1 = [y_1 y_2 y_3 y_4]; 

%F_xy = @(t,r) 3.*exp(-t)-0.4*r; %clc
i = x0;
for i=1:length(x)-1  %funcao der(valor do tempo, vetor com 4 entradas
        k_1 = der2(x(i),y1(:,i));
        k_2 = der2(x(i)+0.5*h,(y1(:,i)+0.5*h*k_1));
        k_3 = der2((x(i)+0.5*h),(y1(:,i)+0.5*h*k_2));
        k_4 = der2((x(i)+h),(y1(:,i)+k_3*h));
        y1(:,i+1) = y1(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h; 
end
plot(x,y1(1,:));grid
title('posição massa 1 por tempo')
figure
plot(x,y1(2,:));grid
title('angulo massa 2 por tempo')
figure
plot(y1(3,:),y1(1,:),y1(2,:),y1(4,:));grid
legend('velocidade x','velocidade alfa')
title('plano de fase velocidade vs posição')
end
