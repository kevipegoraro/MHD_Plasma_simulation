% Funcao dos coeficientes lineares
function c = a(location, state)
global e mi  out ng n0 R0 B0 p0 Vloop % Importar os valores globais usados neste script
gas = 'H2';
tam=length(location.x); % Numero de pontos na mesh
     if tam>1  % Se o calculo é para mais de um ponto, usar todo o vetor interpolado dos campos eletricos e magneticos
         out=campo2(location.x,location.y);
         a12 = first_townsend_coeff(p0,out.Ephi,gas,0);
         Br   = out.Br;
         Bz   = out.Bz;
         Bphi = out.Bphi;
         Ephi = out.Ephi;
         Leff1 = R0*B0./(0.001+(Br.^2 + Bz.^2).^0.5);
     else
            %Caso contrario usar uma aproximacao pela gaussiana para os campos na borda
         a0=0.06;
         di=(0.5/2-location.x);
         s=real(exp(sqrt((di.^2+(location.y).^2))));
         if di==0
             Br=0;
             Bz=0;
         else
             alfa=atan(location.y/di); % Projetando os campos nas direcoes poloidais
             Br   = s*cos(alfa);
             Bz   = s*sin(alfa);
         end
         % Definindo os campos eletricos e magneticos
         Bphi = R0*B0./(R0+location.x*a0);
         Ephi = -Vloop/2/pi./(R0+location.x*a0);
         % Definindo os coeficintes de tranposte
         a12= first_townsend_coeff(p0,Ephi,gas,0);
         Leff1 = R0*B0./(0.001+(Br.^2 + Bz.^2).^0.5);
         
     end
out1 = resistivity_nova(state.u(1,:),state.u(3,:)./state.u(1,:)/e,1);
rho  = out1.eta_par; gas = 'H2';
v_en = 7.89e11*(ng-state.u(1,:))*(1.38e-23)*300;  
v_in = a12.*abs(state.u(2,:))/e./state.u(1,:); 
v_ei = out1.v_ei; 
v_loss = abs(state.u(2,:).*Bphi)/sqrt(Br.^2+ Bphi.^2 +Bz.^2)./state.u(1,:)/e./Leff1; 
%%
c = zeros(4,length(location.x)); % Cria a matri dos coeficientes c com todos valendo zero
c(1,:) = 0;
c(2,:)=v_en+v_in+v_ei-v_loss;
c(3,:)=-2*(e^2)*state.u(1,:).*rho/mi;
c(4,:)=-2*(e^2)*state.u(1,:).*rho/mi;

c=real(c); % Retorna apenas a parte real
