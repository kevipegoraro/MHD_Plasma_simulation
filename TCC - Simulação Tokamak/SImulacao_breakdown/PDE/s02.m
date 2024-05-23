function f = s02(location,state)

N = 2; % Number of equations
nr = length(location.x); % Number of columns
f = zeros(N,nr); % Allocate f

% Now the particular functional form of f
%keyboard
f(1,:) = 3.*real(exp(sqrt(((0.5-location.x).^2+(location.y+0).^2)))); %location.x - location.y + state.u(1,:);
f(2,:) = 5*state.uy(2,:)+location.x.*state.uy(1,:); % + tanh(state.ux(1,:)) + tanh(state.uy(3,:));


%M = length(location.x);
%uinit = zeros(2,M);
%uinit(1,:) = real(exp(sqrt(((0.5-location.x).^2+(location.y+0).^2))));
 