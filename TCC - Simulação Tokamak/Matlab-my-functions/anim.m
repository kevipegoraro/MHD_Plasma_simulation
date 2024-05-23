%Z = peaks;
%figure('Renderer','zbuffer');
%for j = 1:20
 %   surf(sin(2*pi*j/20)*Z,Z)
%    F(j) = getframe;
%end
%movie(F,20,0.1)
function anim(l,x,y1)
figure('Renderer','zbuffer');
for j = 1:length(x)
    %hold on
    plot(y1(1,j),0,'.',y1(1,j)+l*sin(y1(2,j)),-l*cos(y1(2,j)),'.',-l,-1.5*l,'.',1.5*l,1,'.');grid
    %reta da massa1 a massa2:
    %u=[y1(1,j) 0];
    %v=[y1(1,j)+l*sin(y1(2,j)) -l*cos(y1(2,j))];
    %line(u,v)
    F(j) = getframe;
end
movie(F,length(x))
end