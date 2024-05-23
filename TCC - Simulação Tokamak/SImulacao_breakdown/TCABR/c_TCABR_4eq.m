function t = c_TCABR_4eq(location, ~)
D=-0.1;
t = zeros(4,length(location.x));
t(1,:)=D*location.x; %region.y; 