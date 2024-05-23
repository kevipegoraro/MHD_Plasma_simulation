function d=der3(t,x)
% x(1) =>posição y1
% x(2) =>angulo teta
% x(3) =>velocidade dy1/dt
% x(4) =>velocidade angular d(teta)/dt
m1=1; m2=1; c1=1; c2=0.25; L=5; g=9.8; k=20;
I(1,1) = 1; I(1,2)=1; % matriz de inercia
I(2,1) = I(1,2); I(2,2)=1;
C=[0 x(4)*x(1); -2*x(4)*x(1) 0]; %coef de atrito
K=[0 0; 0 0]; %contantes elasticas
F = [L/m1*cos(x(2)); L*sin(x(2))/(m1*x(1)^2)]; %coriolis centrifugo
G = [0 ; 0]; % vetor gravitacional
B = [0;0];
vel=[x(3); x(4)];
pos=[x(1); x(2)];
u = 0; %força externa
v= [(L/m1)*cos(x(2))+x(1)*x(4)^2 ((L/m1)*sin(x(2))-2*x(1)*x(3)*x(4))*(1/x(1)^2)];%inv(I)*(B*u-C*vel-K*pos-F-G);

d(1)=x(3);
d(2)=x(4);
d(3)=v(1);
d(4)=v(2);
d=d';
end