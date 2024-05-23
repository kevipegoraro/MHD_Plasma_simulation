function f = s03(location,state)
%carregando variaveis globais
global e me mi Neq nrpontos out
%definindo a matriz f
tam=length(location.x)
f = zeros(Neq,tam); 
nul = zeros(1,tam);
%Jr(:,:,t) ->state.u(2,:)
%dJr(:,:,t)/dx ->state.ux(2,:)
%n(:,:,t) -> state.u(1,:)
% tempo atual t -> state.time 
% x atual -> location.x
% y atual -> location.y
aux=1:1:tam;
%keyboard
f(1,:) = (state.ux(2,:)+state.uy(3,:)+state.uy(4,:))/e+(e/me)*real(exp(sqrt(((0.5-location.x).^2+(location.y+0).^2)))); %n
f(2,:) = 2*e/me*(out.Bz(aux).*state.u(4,:) - state.u(3,:).*out.Bphi(aux)) - e/me*state.ux(5,:); % f do  J_x r
f(3,:) = e/me*(out.Bphi(aux).*state.u(2,:) - state.u(4,:).*out.Br(aux)) - e/me*state.ux(6,:); % f do  J_y z
f(4,:) = state.u(1,:)*e^2/me.*out.Ephi(aux)+e/me*(out.Br(aux).*state.u(3,:) - state.u(2,:).*out.Bz(aux)); % f do  J_phi
f(5,:) = nul; % f do pe
f(6,:) = nul; % f do pion
f(7,:) = nul; % f do v_en
f(8,:) = nul; % f do v_in
f(9,:) = nul; % f do v_ei
f(10,:) = nul; % f do v_loss
f(11,:) = nul; % f do
f(12,:) = nul; % f do
f(13,:) = nul; % f do
f(14,:) = nul; % f do
f(15,:) = nul; % f do
f(16,:) = nul; % f do
f(17,:) = nul; % f do
f(18,:) = nul; % f do
f(19,:) = nul; % f do


 