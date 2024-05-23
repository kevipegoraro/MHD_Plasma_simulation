function uinit = u0fun(location)
%M = length(location.x);
%uinit = zeros(M);
uinit = exp(sqrt(location.x.^2 + location.y.^2));