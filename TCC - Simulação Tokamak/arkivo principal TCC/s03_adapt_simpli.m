set(0,'defaulttextinterpreter','Latex');
% Funcao do coeficiente f
function f = s03_adapt_simpli(location,state)
% Carregando variaveis globais
global e me mi Neq out rho v_ei ng p0 n0 R0 B0 Vloop 
% Definindo a matriz f
tam=length(location.x);
T=length(out.Br(1,:));
f = zeros(Neq,tam);
for j=1:1:Neq % Checando as derivadas em x das Neq equacoes
    if isnan(state.ux(j,1))
        state.ux(j,:)=0;
    end
    if isnan(state.uy(j,1))
        state.uy(j,:)=0;
    end
end
if state.u(1,:) == 0
    state.u(1,:)=n0;
end
% u(1,:) -> n
% u(2,:) -> J_phi
% u(3,:) -> pe
% u(4,:) -> pion
out1 = resistivity_nova(state.u(1,:),state.u(3,:)./state.u(1,:)/e,1);
rho  = out1.eta_par; gas = 'H2';

     if tam>1 % Se o calculo e para mais de um ponto, usar todo o vetor interpolado dos campos eletricos e magneticos
         out=campo2(location.x,location.y);
         a12 = first_townsend_coeff(p0,out.Ephi,gas,0);
         Br   = out.Br;
         Bz   = out.Bz;
         Bphi = out.Bphi;
         Ephi = out.Ephi;
         Leff1 = R0*B0./(0.001+(Br.^2 + Bz.^2).^0.5);
     else
            % Caso contrario usar uma aproximacao pela gaussiana para os campos na borda
         a0=0.06;
         di=(0.5/2-location.x);
         s=real(exp(sqrt((di.^2+(location.y).^2))));
         if di==0
             Br=0;
             Bz=0;
         else
             alfa=atan(location.y/di);
             Br   = s*cos(alfa);
             Bz   = s*sin(alfa);
         end
         Bphi = R0*B0./(R0+location.x*a0);
         Ephi = -Vloop/2/pi./(R0+location.x*a0);
         a12= first_townsend_coeff(p0,Ephi,gas,0);
         Leff1 = R0*B0./(0.001+(Br.^2 + Bz.^2).^0.5);
     end
v_en = 7.89e11*(ng-state.u(1,:))*(1.38e-23)*300;  
v_in = a12.*abs(state.u(2,:))/e./state.u(1,:); 
v_ei = out1.v_ei; % f do  
v_loss = abs(state.u(2,:).*Bphi)/sqrt(Br.^2+Bphi.^2 + Bz.^2)./state.u(1,:)/e./Leff1; 

% Densidade de particulas 
f(1,:) = (state.ux(2,:)+state.uy(2,:))/e+v_in-v_loss;  %f do n 

% Densidade de corrente J_phi
f(2,:) = state.u(1,:)*e^2/me.*Ephi; % f do  J_phi

% Pressao pe e pion
f(3,:) = 3/2*(1 + (2*v_en+v_in-v_loss)./v_ei/2).*rho.*(state.u(2,:).^2)-2*state.u(1,:)*e^2/mi.*rho.*state.u(4,:);    %f do pe

f(4,:) = 2*state.u(1,:)*e^2/mi.*rho.*state.u(3,:); % f do pion

f=real(f); % Retorna a parte real de f
 