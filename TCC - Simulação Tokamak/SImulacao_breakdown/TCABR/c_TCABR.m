function t = c_TCABR(location, ~)
D=0.4;
t = zeros(6,length(location.x));
t(1,:)=D*location.x; %region.y; 