%função condições iniciais
function uinit = inifunc_simpli(location, state)
%primeira edp U= zeros(Neq,1); U(1)=n0; U(5)=e*Te0*n0; U(6)=U(5);
global Neq Te0 n0 e out1
%keyboard
uinit = zeros(Neq,length(location.x)); %cria com todos zeros
uinit(1,:) = n0;
uinit(2,:)=-1e-6;
uinit(3,:)=e*Te0*n0;
uinit(4,:)=e*Te0*n0;
uinit(5,:)= 0;
uinit(6,:)= out1.v_ei;
uinit(7,:)= 0;

%
% u(1,:) -> n        state.u(1,:)
% u(2,:) -> J_phi    state.u(2,:)
% u(3,:) -> pe       state.u(3,:)
% u(4,:) -> pion     state.u(4,:)
% u(5,:) -> v_in     state.u(5,:)
% u(6,:) -> v_ei     state.u(6,:)
% u(7,:) -> v_loss   state.u(7,:)