function uinit = u0fun(location)
%primeira edp
M = length(location.x);
uinit = ones(2,M);
uinit(1,:) = 10; %4 + location.x.^2 + location.y.^2;