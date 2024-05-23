function f = s03(location,state)
%carregando variaveis globais
global e me mi G 
G.location.x=location.x;
G.location.y=location.y;
G.tam=length(G.location.x);

f = zeros(4,G.tam);

for j=1:1:4 %checando as derivadas em x das Neq equações
    if isnan(state.ux(j,1))
        state.ux(j,:)=0;
    end
    if isnan(state.uy(j,1))
        state.uy(j,:)=0;
    end
end
%ind = find(state.u(1,:)<G.n0); state.u(1,ind) = G.n0;
%ind = find(state.u(3,:)<e*G.Te0*G.n0); state.u(3,ind) =e*G.Te0*G.n0;
%ind = find(state.u(4,:)<e*G.Te0*G.n0); state.u(4,ind) =e*G.Te0*G.n0;

% u(1,:) -> n
% u(2,:) -> J_phi, J_r=0
% u(3,:) -> pe
% u(4,:) -> pion

G.out1 = resistivity_nova(state.u(1,:),state.u(3,:)./state.u(1,:)/e,1);
G.rho  = G.out1.eta_par; 

out=RegulaCampo; % Ajuste dos campos eletromagneticos para o numero de pontos contidos em location 
if G.ligacorrente % Se ligacorrente for 1, inclui-se o campo magnetico gerado pela corrente de plasma
   % corrente = 0; area = 10000000*3.1415*0.06^2; % Calculando a corrente total
   % for ll=1:length(state.u(1,:)) corrente = corrente+state.u(2,ll); end
   corrente =  10e5*3.1415*0.06^2;
   out.Bplx = out.Bplx*corrente;
    out.Bply = out.Bply*corrente;
    out.Br = out.Br+out.Bplx;
    out.Bz = out.Bz+out.Bply;
end
%keyboard
v_en = G.constante*(G.ng-state.u(1,:));  
v_in = out.a1.*abs(state.u(2,:))/e./state.u(1,:); 
v_ei = G.out1.v_ei;   
v_loss = abs(state.u(2,:).*out.Bphi)/sqrt(out.Br.^2+out.Bphi.^2 + out.Bz.^2)./state.u(1,:)/e./out.Leff; 
s=G.ligas*10e10*exp(-1000*((-G.R0+location.x).^2+(location.y).^2)); %Aproximacao por uma Gaussiana
v_in= G.ligacoef*v_in*10e7; % Aumento intensional da quantidade de ionizacao e perda
v_loss= G.ligacoef*v_loss*10e7;
if G.plotton
if G.tam>1  % Se plotton for 1 e salvo uma substrutura com os valores de v_in-v_loss para plotagem posterior 
        eval(['G.plottv' num2str(G.cplott) '=v_in;'])
        eval(['G.plottl' num2str(G.cplott) '=v_loss;'])
       % if G.t != state.time
        eval(['G.plottcampoplasma' num2str(G.cplott) '=sqrt(out.Br.^2+out.Bz.^2);'])
        %end
       % G.t=state.time;
        G.cplott=G.cplott+1;
end
end
%
%densidade de particulas 

au = (v_in-v_loss+s); % 
%au = s;
f(1,:) = location.x.*au;
%densidade de corrente J_phi
f(2,:) = state.u(1,:)*e^2/me.*out.Ephi; % f do  J_phi
%pressão pe e pion
f(3,:) = 3/2*(1 + (2*v_en+au)./v_ei/2).*G.rho.*(state.u(2,:).^2)+2*state.u(1,:)*e^2/mi.*G.rho.*state.u(4,:);    %f do pe
f(4,:) = 2*state.u(1,:)*e^2/mi.*G.rho.*state.u(3,:); % f do pion

f=real(f);
 