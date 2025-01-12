function c = a01_teste(location, state)
global e mi G
out1 = resistivity_nova(state.u(1,:),state.u(3,:)./state.u(1,:)/e,1);
rho  = out1.eta_par;
s=state.time*10e12*exp(-1000*((-G.R0+location.x).^2+(location.y).^2));
%s=G.ligas*state.time*10e10*exp(-1000*((-G.R0+location.x).^2+(location.y+0.04-state.time*3).^2))+...
 %   (0.02-state.time)*10e10*exp(-1000*((-G.R0+location.x).^2+(location.y-0.04+state.time*3).^2))+...
  %  (0.02-state.time)*10e10*exp(-1000*((-G.R0+location.x+0.04-state.time*3).^2+(location.y).^2))+...
   % state.time*10e10*exp(-1000*((-G.R0+location.x-0.04+state.time*3).^2+(location.y).^2));
c = zeros(4,length(location.x)); % Cria uma matriz com todos os pontos sendo zero
c(1,:) = 0;
c(2,:)= s; % Calculo pelos coeficiente de ionizacao e de perda if G.tam>1 keyboard; end
c(3,:)=2*(e^2)*state.u(1,:).*rho/mi;
c(4,:)=2*(e^2)*state.u(1,:).*rho/mi;
c=real(c);
