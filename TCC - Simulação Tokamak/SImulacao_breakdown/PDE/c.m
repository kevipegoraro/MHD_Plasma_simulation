function t = c(location, state)
%keyborad
D=-0.02;
t = zeros(4,length(location.x));
t(1,:)=D*location.x; 
end