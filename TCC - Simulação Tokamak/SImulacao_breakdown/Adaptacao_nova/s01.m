function f =s01(location,state)
nr = length(location.x); % Number of columns
f = zeros(nr); % Allocate f
% Now the particular functional form of f
f(1,:) = real(exp((0.03+location.x).^2-location.y.^2)/2)) + state.u(1,:);
%f(2,:) = 1 + tanh(state.ux(1,:)) + tanh(state.uy(3,:));
%f(3,:) = (5 + state.u(3,:)).*sqrt(location.x.^2 + location.y.^2);