% Funcao condicoes iniciais
function uinit = inifuncadapt(location, state)
global Neq Te0 n0 e 
r=1-(location.x.^2+location.y.^2).^2; % Garante que J_phi, Pe, Pi na borda sera 0.
uinit = zeros(Neq,length(location.x)); % Cria com todos zeros
uinit(1,:) = n0; % A condicao inicial da densidade de particulas em todos os pontos sera n0
uinit(2,:)=-r*1e-6;
uinit(3,:)=r*e*Te0*n0;
uinit(4,:)=r*e*Te0*n0;
