function uinit = inifuncadapt01(location, state)
global e me G
Ephi =  -G.Vloop/2/pi./(G.R0+location.x*G.a0); 
ta=e^2/G.somatranpost/me;
uinit = zeros(4,length(location.x)); %cria com todos zeros
uinit(1,:) = G.n0; %location.x.*
uinit(2,:)=Ephi.*ta.*uinit(1,:);
uinit(3,:)=e*G.Te0*uinit(1,:);
uinit(4,:)=e*G.Te0*uinit(1,:);
