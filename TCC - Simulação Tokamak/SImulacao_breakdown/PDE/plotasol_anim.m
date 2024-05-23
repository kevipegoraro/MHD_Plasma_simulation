function plotasol_anim(model,sol)
global G
s=length(sol(1,1,:));
figure(6); 
a1=max(sol(:,1,:)); 
a2=max(-sol(:,2,:)); 
a3=max(sol(:,3,:));
a4=max(sol(:,4,:)); 
e= 1.6e-19;
for kk=3:7:s
subplot(2,2,1)
pdeplot(model,'XYData',sol(:,1,kk))
X2=['n(dt*' num2str(kk) ')'];
title(X2)
colormap(jet);
caxis([G.n0 a1(kk)])
axis equal 

subplot(2,2,2)
pdeplot(model,'XYData',sol(:,2,kk))
X2=['J_p_h_i(dt*' num2str(kk) ')'];
title(X2)
colormap(jet);
caxis([0 a2(kk)])
axis equal 

subplot(2,2,3)
pdeplot(model,'XYData',sol(:,3,kk))
X2=['pe(dt*' num2str(kk) ')'];
title(X2)
colormap(jet);
caxis([e*G.Te0*G.n0 a3(kk)])
axis equal 

subplot(2,2,4)
pdeplot(model,'XYData',sol(:,4,kk))
X2=['pi(dt*' num2str(kk) ')'];
title(X2)
colormap(jet);
caxis([e*G.Te0*G.n0 a4(kk)])
axis equal

pause(0.3);
end
end