function f = s03_adapt(location,state)
%carregando variaveis globais
global e me mi Neq out constante ng rho v_ei 
%definindo a matriz f
tam=length(location.x);
f = zeros(Neq,tam);
aux=1:1:tam;
for j=1:1:Neq %checando as derivadas em x das Neq equações
    if isnan(state.ux(j,1))
        state.ux(j,:)=0;
    end
    if isnan(state.uy(j,1))
        state.uy(j,:)=0;
    end
end
% u(1,:) -> n
% u(2,:) -> J_phi
% u(3,:) -> pe
% u(4,:) -> pion

%coeficientes de transporte
%out1 = resistivity_nova(state.u(1,:),state.u(3,:)./state.u(1,:)/e,1);
%rhoP=rho;
%rho  = out1.eta_par; 
%if isnan(rho(1)) rho=rhoP; end

v_en = constante*(ng-state.u(1,:));
%v_ei=out1.v_ei
%ind = find(v_ei<0); v_ei(ind) = 0.001;
%if isnan(v_ei) v_ei=0.001; end

%v_ion = a1(aux).*abs(state.u(2,:))/e./state.u(1,:); ind = find(v_ion<0); v_ion(ind) = 0.001;
%auxx=sqrt(out.Br(aux).^2 + out.Bphi(aux).^2 + out.Bz(aux).^2);
% ind = find(auxx<0); uuxx(ind) = 0.001;
%v_loss = abs(state.u(2,:).*out.Bphi(aux))/auxx./state.u(1,:)/e./Leff(aux); 
s=real(exp(sqrt(((0.5/2-location.x).^2+(location.y+0).^2))));

%densidade de particulas  dn/dt - div(J / e) = Se/me + D nabla^2( n ).
f(1,:) = (state.ux(2,:)+state.uy(2,:))/e+s*e/me; %v_ion-v_loss; %f do n 
%ind = find(f(1,:)<n0); f(1,ind) = n0;
%densidade de corrente J_phi
f(2,:) = state.u(1,:)*e^2/me.*out.Ephi(aux); % f do  J_phi

%pressão pe e pion
f(3,:) = 3/2*(1 + (2*v_en(aux)+s)./v_ei/2).*rho.*(state.u(2,:).^2)-2*state.u(1,:)*e^2/mi.*rho.*state.u(4,:);    %f do pe

f(4,:) = 2*state.u(1,:)*e^2/mi.*rho.*state.u(3,:); % f do pion

f=real(f);
 