function f = s03_Gaussiana(location,state)
%carregando variaveis globais
global e me mi G 
G.tam=length(G.location.x);
if state.u(1,:) == 0 state.u(1,:) = G.n0; end
f = zeros(4,G.tam);
for j=1:1:4 %checando as derivadas em x das Neq equações
    if isnan(state.ux(j,1))
        state.ux(j,:)=0;
    end
    if isnan(state.uy(j,1))
        state.uy(j,:)=0;
    end
end

%keyboard
% u(1,:) -> n
% u(2,:) -> J_phi, J_r=0
% u(3,:) -> pe
% u(4,:) -> pion
G.out1 = resistivity_nova(state.u(1,:),state.u(3,:)./state.u(1,:)/e,1);
G.rho  = G.out1.eta_par; 
v_en = G.constante*(G.ng-state.u(1,:));  
v_ei = G.out1.v_ei;    
s=10e8*exp(-1000*((location.x+G.a0/4).^2+(location.y).^2));
Ephi = -G.Vloop/2/pi./location.x;

%densidade de particulas 
f(1,:) = location.x.*s;
%densidade de corrente J_phi
f(2,:) = state.u(1,:)*e^2/me.*Ephi; % f do  J_phi
%pressão pe e pion
f(3,:) = 3/2*(1 + (2*v_en+s)./v_ei/2).*G.rho.*(state.u(2,:).^2)+2*state.u(1,:)*e^2/mi.*G.rho.*state.u(4,:); %f do pe
f(4,:) = 2*state.u(1,:)*e^2/mi.*G.rho.*state.u(3,:); % f do pion
f
end