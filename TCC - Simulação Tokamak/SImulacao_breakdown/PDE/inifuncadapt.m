function uinit = inifuncadapt(location, state)
global e me G
Ephi =  -(G.Vloop/2/pi)./location.x; 
ta=e^2/G.somatranpost/me;
uinit = zeros(4,length(location.x)); %cria com todos zeros
uinit(1,:) = G.n0; 
uinit(2,:)=ta*G.n0.*Ephi;
uinit(3,:)=e*G.Te0*uinit(1,:);
uinit(4,:)=e*G.Te0*uinit(1,:);
end