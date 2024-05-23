function c = a_TCABR_4eq(location, state)
global e mi G
s=G.s1*exp(-G.s2*((G.s3+location.x).^2+(location.y).^2)); %state.time
v_en= G.constante*(G.ng-state.u(1,:)); 
if state.u(1,:) == 0 aux=state.u(3,:)./G.n0/e; else aux=state.u(3,:)./state.u(1,:)/e; end
out1 = resistivity_TCABR(state.u(1,:),aux,1);
rho  = out1.eta_par; 
c = zeros(4,length(location.x)); % Cria uma matriz com todos os pontos sendo zero
c(2,:)=v_en+s-e^2*rho.*state.u(1,:);
c(3,:)=2*(e^2)*state.u(1,:).*rho/mi; 
c(4,:)=2*(e^2)*state.u(1,:).*rho/mi;
c=real(c);
