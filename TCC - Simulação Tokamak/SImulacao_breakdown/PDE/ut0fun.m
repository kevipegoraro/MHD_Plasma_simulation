function utinit = ut0fun(location)
%segunda edp
M = length(location.x);
utinit = ones(2,M);
utinit(2,:) = 1;%sin(location.x.*location.y);
