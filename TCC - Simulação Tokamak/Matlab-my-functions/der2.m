function d=der2(t,x)
% x(1) =>posição y1
% x(2) =>angulo teta
% x(3) =>velocidade dy1/dt
% x(4) =>velocidade angular d(teta)/dt
m1=2; m2=1; c1=1; c2=0.25; L=5; g=9.8; k=20;
I(1,1) = m1+m2; I(1,2)=m2*L*cos(x(2)); % matriz de inercia
I(2,1) = I(1,2); I(2,2)=m2*L^2;
C=[c1 0; 0 c2]; %coef de atrito
K=[k 0; 0 0]; %contantes elasticas
F = [-m2*L*sin(x(2))*x(4)^2; 0]; %coriolis centrifugo
G = [0 ; m2*g*L*sin(x(2))]; % vetor gravitacional
B = [1;0];
vel=[x(3); x(4)];
pos=[x(1); x(2)];
u = 0; %força externa
v=inv(I)*(B*u-C*vel-K*pos-F-G);
d(1)=x(3);
d(2)=x(4);
d(3)=v(1);
d(4)=v(2);
d=d';
end