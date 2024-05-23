function t = c01(location, state)
D=0.03;
t = zeros(4,length(location.x));
t(1,:)=D*location.x; %region.y; 