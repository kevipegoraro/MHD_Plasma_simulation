function anim_q(l,x,y1)
figure('Renderer','zbuffer');
for j = 1:length(x)
    %hold on
    plot(y1(1,j)*cos(y1(2,j)),y1(1,j)*sin(y1(2,j)),'.',-l,l,'.',l,-l,'.');grid
    %reta da massa1 a massa2:
    %u=[y1(1,j) 0];
    %v=[y1(1,j)+l*sin(y1(2,j)) -l*cos(y1(2,j))];
    %line(u,v)
    F(j) = getframe;
end
movie(F,length(x))
end