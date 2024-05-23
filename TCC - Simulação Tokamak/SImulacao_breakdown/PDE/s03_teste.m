function f = s03_teste(location,state)
global e me mi G 
G.tam=length(location.x);
f = zeros(4,G.tam);
for j=1:1:4 
    if isnan(state.ux(j,1))
        state.ux(j,:)=0;
    end
    if isnan(state.uy(j,1))
        state.uy(j,:)=0;
    end
end 
s=state.time*10e12*exp(-1000*((-G.R0+location.x).^2+(location.y).^2));
%s=state.time*10e10*exp(-1000*((-G.R0+location.x).^2+(location.y+0.04-state.time*3).^2))+...
  %  (0.02-state.time)*10e10*exp(-1000*((-G.R0+location.x).^2+(location.y-0.04+state.time*3).^2))+...
  %  (0.02-state.time)*10e10*exp(-1000*((-G.R0+location.x+0.04-state.time*3).^2+(location.y).^2))+...
   % state.time*10e10*exp(-1000*((-G.R0+location.x-0.04+state.time*3).^2+(location.y).^2));
%keyboard
    Ephi = -G.Vloop/2/pi./(G.R0+location.x*G.a0);
v_en= G.constante*(G.ng-state.u(1,:)); 
out1 = resistivity_nova(state.u(1,:),state.u(3,:)./state.u(1,:)/e,1);
rho  = out1.eta_par;
%v_ei=out1.v_ei;
f(1,:) = location.x.*s;
f(2,:) = state.u(1,:)*e^2/me.*Ephi; 
f(3,:) = 2*state.u(1,:)*e^2/mi.*rho.*state.u(4,:)+v_en;
%max(f(3,:))  f(3,1) f(3,end)
%G.n0*e*G.Te0; %3/2*(1 + (2*v_en+s)./v_ei/2).*G.rho.*(state.u(2,:).^2)+2*state.u(1,:)*e^2/mi.*G.rho.*state.u(4,:);
f(4,:) = 2*state.u(1,:)*e^2/mi.*G.rho.*state.u(3,:); 
f=real(f);
 