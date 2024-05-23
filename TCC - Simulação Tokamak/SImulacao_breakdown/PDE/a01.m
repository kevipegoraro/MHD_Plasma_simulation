function c = a01(location, state)
global e mi G

G.location.x=location.x;
G.location.y=location.y;
G.tam=length(G.location.x);

out=RegulaCampo;
if G.ligacorrente % Se ligacorrente for 1, inclui-se o campo magnetico gerado pela corrente de plasma
   % corrente = 0; area = 10000000*3.1415*0.06^2; % Calculando a corrente total
   % for ll=1:length(state.u(1,:)) corrente = corrente+state.u(2,ll); end
   corrente =  10e5*3.1415*0.06^2;
   out.Bplx = out.Bplx*corrente;
    out.Bply = out.Bply*corrente;
    out.Br = out.Br+out.Bplx;
    out.Bz = out.Bz+out.Bply;
end

G.out1 = resistivity_nova(state.u(1,:),state.u(3,:)./state.u(1,:)/e,1);
rho  = G.out1.eta_par;

v_en = G.constante*(G.ng-state.u(1,:)); 

v_in = out.a1.*abs(state.u(2,:))/e./state.u(1,:);  
v_ei = G.out1.v_ei; 
v_loss = abs(state.u(2,:).*out.Bphi)/sqrt(out.Br.^2+out.Bphi.^2 + out.Bz.^2)./state.u(1,:)/e./out.Leff; 
s=G.ligas*10e10*exp(-1000*((-G.R0+location.x).^2+(location.y).^2)); %Aproximacao por uma Gaussiana
v_in= G.ligacoef*v_in*10e7; % Aumento intensional da quantidade de ionizacao e perda
v_loss= G.ligacoef*v_loss*10e7;

c = zeros(4,length(location.x)); % Cria uma matriz com todos os pontos sendo zero
c(1,:) = 0;
c(2,:)=v_en+v_in+v_ei-v_loss+s; % Calculo pelos coeficiente de ionizacao e de perda if G.tam>1 keyboard; end
c(3,:)=2*(e^2)*state.u(1,:).*rho/mi;
c(4,:)=2*(e^2)*state.u(1,:).*rho/mi;

c=real(c);
