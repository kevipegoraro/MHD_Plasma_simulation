function f = f_TCABR_4eq(location,state)
global e me mi G 
G.tam=length(location.x);
G.location.x=location.x;
G.location.y=location.y;
f = zeros(4,G.tam);
if state.u(1,:)==0 state.u(1,:)=G.n0; end
for j=1:1:4 
    if isnan(state.ux(j,1))
        state.ux(j,:)=0;
    end
    if isnan(state.uy(j,1))
        state.uy(j,:)=0;
    end
end

%out=RegulaCampo; % Ajuste dos campos eletromagneticos para o numero de pontos contidos em location 
%if G.ligacorrente % Se ligacorrente for 1, inclui-se o campo magnetico gerado pela corrente de plasma
   % corrente = 0; area = 10000000*3.1415*0.06^2; % Calculando a corrente total
   % for ll=1:length(state.u(1,:)) corrente = corrente+state.u(2,ll); end
  %  corrente =  10e5*3.1415*0.06^2;
 %   out.Bplx = out.Bplx*corrente;
%    out.Bply = out.Bply*corrente;
  %  out.Br = out.Br+out.Bplx;
 %   out.Bz = out.Bz+out.Bply;
%end
%Bx = out.Br; Bz = out.Bz;
%keyboard
s=G.s1*exp(-G.s2*((G.s3+location.x).^2+(location.y).^2));
Ephi = -G.Vloop/2/pi./(G.R0+location.x); 
if state.u(1,:) == 0 aux=state.u(3,:)./G.n0/e; else aux=state.u(3,:)./state.u(1,:)/e; end
out1 = resistivity_TCABR(state.u(1,:),aux,1);
rho  = out1.eta_par;
v_en = G.constante*(G.ng-state.u(1,:));
%Bz.*state.u(3,:).*state.u(1,:)+state.ux(5,:)/me^2;
%Bx.*state.u(3,:).*state.u(1,:)+state.uy(5,:)/me^2;
f(1,:) = location.x.*s;
f(2,:) = -e.*Ephi+0*(state.ux(3,:)+state.uy(3,:))./state.u(1,:);  
f(3,:) = 2*e^2/mi*state.u(1,:).*rho.*state.u(4,:) ...
    +2/3*state.u(1,:).*(state.ux(3,:)+state.uy(3,:)) ...
    +(state.u(2,:).^2).*state.u(1,:).*(1/2*me*(2*v_en+s)-e^2*rho.*state.u(1,:));
%if G.tam>1 
%    keyboard; 
%max(f(3,:))  
%f(3,1) 
%f(3,end)
%end
f(4,:) = 2*state.u(1,:)*e^2/mi.*rho.*state.u(3,:); 
f=real(f);
 