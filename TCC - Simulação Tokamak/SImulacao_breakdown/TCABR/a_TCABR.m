function c = a_TCABR(location, state)
global e mi G
s=G.s1*exp(-G.s2*((G.s3+location.x).^2+(location.y).^2)); %state.time
v_en= G.constante*(G.ng-state.u(1,:)); 
out1 = resistivity_nova(state.u(1,:),state.u(3,:)./state.u(1,:)/e,1);
rho  = out1.eta_par; 
c = zeros(6,length(location.x)); % Cria uma matriz com todos os pontos sendo zero
c(3,:)=v_en+s-e^2*rho.*state.u(1,:);
c(5,:)=2*(e^2)*state.u(1,:).*rho/mi; 
c(6,:)=2*(e^2)*state.u(1,:).*rho/mi;
% ----teste----
c(2,:)=G.uxcoeficientes;
c(4,:)=G.uxcoeficientes;
% ----teste----
c=real(c);
