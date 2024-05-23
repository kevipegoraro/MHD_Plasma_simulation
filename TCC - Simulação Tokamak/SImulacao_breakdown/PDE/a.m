function c = a(location, state)
global e mi G
G.tam=length(location.x);
%ind = find(state.u(1,:)<G.n0); state.u(1,ind) = G.n0;
%ind = find(state.u(3,:)<e*G.Te0*G.n0); state.u(3,ind) =e*G.Te0*G.n0;
%ind = find(state.u(4,:)<e*G.Te0*G.n0); state.u(4,ind) =e*G.Te0*G.n0;
out1 = resistivity_nova(state.u(1,:),state.u(3,:)./state.u(1,:)/e,1);
v_en = G.constante*(G.ng-state.u(1,:));  
if state.time > 25*G.dt s=0; else s=10e11/(G.dt*25)*state.time*exp(-10000*((-G.R0+location.x+G.a0/4).^2+(location.y).^2)); end
%s=10e11*exp(-10000*((-G.R0+location.x+G.a0/4).^2+location.y.^2));
c = zeros(4,length(location.x)); 
K=2*(e^2)*state.u(1,:).*out1.eta_par/mi;
c(1,:) = 0;
c(2,:)=v_en+out1.v_ei+s;
c(3,:)=K;
c(4,:)=K;
c=real(c);
end