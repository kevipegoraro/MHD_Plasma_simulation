function f= s03_Gaussiana(location, state)
global e me mi G 
G.tam=length(location.x);
%ind = find(state.u(1,:)<G.n0); state.u(1,ind) = G.n0;
%ind = find(state.u(3,:)<e*G.Te0*G.n0); state.u(3,ind) =e*G.Te0*G.n0;
%ind = find(state.u(4,:)<e*G.Te0*G.n0); state.u(4,ind) =e*G.Te0*G.n0;
if state.time > 25*G.dt s=0; else s=10e11/(G.dt*25)*state.time*exp(-10000*((-G.R0+location.x+G.a0/4).^2+(location.y).^2)); end %Aproximacao por uma Gaussiana
f = zeros(4,G.tam);
out1 = resistivity_nova(state.u(1,:),state.u(3,:)./state.u(1,:)/e,1);
rho  = out1.eta_par; 
Ephi =  -(G.Vloop/2/pi)./location.x; 
%keyboard
v_en = G.constante*(G.ng-state.u(1,:));  
v_ei = out1.v_ei;   
if G.tam>1
    if G.plotton 
        G.plott1=s; G.plotton=0; 
    end
end
%densidade de particulas 
f(1,:) = location.x.*s;
%densidade de corrente J_phi
f(2,:) = state.u(1,:)*e^2/me.*Ephi; % f do  J_phi
%pressão pe e pion
f(3,:) = 3/2*(1 + (2*v_en+s)./v_ei/2).*rho.*(state.u(2,:).^2)+2*state.u(1,:)*e^2/mi.*rho.*state.u(4,:);    %f do pe
f(4,:) = (2*state.u(1,:)*e^2/mi).*rho.*state.u(3,:); % f do pion

f=real(f);
 