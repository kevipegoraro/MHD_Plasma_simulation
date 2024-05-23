function d=der3(t,x)
% x(1) =>posição x
% x(2) =>posição y
% x(3) =>velocidade dx/dt
% x(4) =>velocidade dy/dt
m1=1; f1=4; f2=1; ii=1; g=9.8; 
b = 2 ; %força b(t) externa
v = [b*sin(2*ii*(x(1)+f2))/m1 b*cos(2*ii*(x(1)+f2))/m1+g]; 

d(1)=x(3);
d(2)=x(4);
d(3)=v(1);
d(4)=v(2);
d=d';
end