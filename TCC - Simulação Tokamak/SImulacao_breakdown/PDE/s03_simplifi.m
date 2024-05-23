function f = s03(location,state)
%carregando variaveis globais
global e me mi Neq out Leff ng a1 n0 constante auxx
%definindo a matriz f
tam=length(location.x);
f = zeros(Neq,tam);
aux=1:1:tam;
%Info = [state.time max(state.ux(1,:)) max(state.ux(2,:)) max(state.ux(3,:)) max(state.ux(4,:)) max(state.ux(5,:)) max(state.ux(6,:))]

if state.u(1,1) == 0
    state.u(1,:)=n0;
end
out1 = resistivity_nova(state.u(1,:),state.u(3,:)./state.u(1,:)/e,1);
rho  = out1.eta_par;
%keyboard

% tempo atual t -> state.time 
% x atual -> location.x
% y atual -> location.y

%
% u(1,:) -> n        state.u(1,:)
% u(2,:) -> J_phi    state.u(2,:)
% u(3,:) -> pe       state.u(3,:)
% u(4,:) -> pion     state.u(4,:)
% u(5,:) -> v_in     state.u(5,:)
% u(6,:) -> v_ei     state.u(6,:)
% u(7,:) -> v_loss   state.u(7,:)
%tirei o maldito v_en
v_en = constante*(ng-state.u(1,:));
%keyboard
%densidade de particulas
f(1,:) = (state.ux(2,:)+state.uy(2,:))/e+state.u(5,:)-state.u(7,:); %f do n 
ind = find(f(1,:)<n0);
f(1,ind) = n0;
%densidade de corrente
f(2,:) = state.u(1,:)*e^2/me.*out.Ephi(aux);

%pressoes 
f(3,:) =  3/2*(1 + (2*v_en+state.u(5,:)-state.u(7,:))./(state.u(6,:)+0.01)/2).*rho.*(state.u(2,:).^2);  %-2*state.u(1,:)*e^2/mi.*rho.*state.u(4,:);
f(4,:) =  2*state.u(1,:)*e^2/mi.*rho.*state.u(3,:);  
%coeficientes de transporte                                      
f(5,:) = a1(aux).*abs(state.u(2,:))/e./state.u(1,:); % f do v_in
f(6,:) = out1.v_ei; % f do v_ei 
f(7,:) = abs(state.u(2,:).*out.Bphi(aux))/auxx./state.u(1,:)/e./Leff(aux);  % f do v_loss
 %size(f)  
 f=real(f);
 