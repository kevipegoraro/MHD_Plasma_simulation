function c = a_TCABR(location, state)
global e mi 
out1 = resistivity_nova(state.u(1,:),state.u(3,:)./state.u(1,:)/e,1);
rho  = out1.eta_par;
s=state.time*10e10*exp(-100*((-0.65+location.x).^2+(location.y).^2));
c = zeros(4,length(location.x)); % Cria uma matriz com todos os pontos sendo zero
c(1,:) = 0;
c(2,:)= s; % Calculo pelos coeficiente de ionizacao e de perda if G.tam>1 keyboard; end
c(3,:)=2*(e^2)*state.u(1,:).*rho/mi;
c(4,:)=2*(e^2)*state.u(1,:).*rho/mi;
c=real(c);
