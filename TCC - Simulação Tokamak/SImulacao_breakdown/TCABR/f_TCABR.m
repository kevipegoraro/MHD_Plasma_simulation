function f = f_TCABR(location,state)
global e me mi G 
G.tam=length(location.x);
G.location.x=location.x;
G.location.y=location.y;
f = zeros(4,G.tam);

for j=1:1:4 
    if isnan(state.ux(j,1))
        state.ux(j,:)=0;
    end
    if isnan(state.uy(j,1))
        state.uy(j,:)=0;
    end
end

out=RegulaCampo; % Ajuste dos campos eletromagneticos para o numero de pontos contidos em location 
if G.ligacorrente 
    corrente =  10e5*3.1415*0.06^2;
    out.Bplx = out.Bplx*corrente;
    out.Bply = out.Bply*corrente;
    out.Br = out.Br+out.Bplx;
    out.Bz = out.Bz+out.Bply;
end
Bx = out.Br;
Bz = out.Bz;

s=G.s1*exp(-G.s2*((G.s3+location.x).^2+(location.y).^2));
Ephi = -G.Vloop/2/pi./(G.R0+location.x*G.a0); 
out1 = resistivity_nova(state.u(1,:),state.u(3,:)./state.u(1,:)/e,1);
rho  = out1.eta_par;
f(1,:) = location.x.*s;
f(2,:) = Bz.*state.u(3,:).*state.u(1,:)+state.ux(5,:)/me^2;
f(3,:) = -e.*Ephi-(state.ux(3,:)+state.uy(3,:))./state.u(1,:);
f(4,:) = Bx.*state.u(3,:).*state.u(1,:)+state.uy(5,:)/me^2;
f(5,:) = 2*state.u(1,:)*e^2/mi.*rho.*state.u(4,:);
%max(f(3,:))  f(3,1) f(3,end)
%G.n0*e*G.Te0; %+2*state.u(1,:)*e^2/mi.*G.rho.*state.u(4,:);
f(6,:) = 2*state.u(1,:)*e^2/mi.*G.rho.*state.u(3,:); 
f=real(f);
 