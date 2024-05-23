function uinit = ini_TCABR(location, ~)
global e me G
Ephi = -G.Vloop/2/pi./(G.R0+location.x*G.a0); 
ta=e/G.somatranpost/me;
uinit = zeros(6,length(location.x)); %cria com todos zeros
uinit(1,:) = G.n0; %location.x.*
% ----teste----
%keyboard
uinit(2,:) = G.uxcoeficientes2;
uinit(4,:) = G.uxcoeficientes2;
% ----teste----
uinit(3,:)=Ephi.*ta;
uinit(5,:)=e*G.Te0*uinit(1,:);
uinit(6,:)=e*G.Te0*uinit(1,:);
