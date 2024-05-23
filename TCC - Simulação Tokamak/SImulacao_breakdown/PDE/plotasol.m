function plotasol(model,sol,c)
global G e
s=length(sol(1,1,:));
T=[3 round(s/4,0) round(3*s/4,0) s-1];
for i=c
    if i==1
        X=['n(:,'];
        amax=G.n0*1.001;%max(max(sol(:,i,:)))/5
        amim=G.n0
    elseif i==2
        X=['J_p_h_i(:,'];
        amax=5e-9/5;%max(max(-sol(:,i,:)));
        amim=-max(max(sol(:,i,:)));
    elseif i==3
        X=['Pe(:,'];
        amax=e*G.Te0*G.n0*1.01; %max(max(sol(:,i,:)))
        amim=e*G.Te0*G.n0
    elseif i==4
        X=['Pi(:,'];
        amax=e*G.Te0*G.n0*1.0001
        amim=e*G.Te0*G.n0
    end
figure(i)
subplot(2,2,1)
pdeplot(model,'XYData',sol(:,i,T(1)))
X2=[X 'dt*' num2str(T(1)) ')'];
title(X2)
colormap(jet);% caxis([amim amax]);
axis equal 
subplot(2,2,2)
pdeplot(model,'XYData',sol(:,i,T(2)))
X2=[X 'dt*' num2str(T(2)) ')'];
title(X2)
colormap(jet); %caxis([amim amax]);
axis equal 
subplot(2,2,3)
pdeplot(model,'XYData',sol(:,i,T(3)))
X2=[X 'dt*' num2str(T(3)) ')'];
title(X2)
colormap(jet); %caxis([amim amax]);
axis equal 
subplot(2,2,4)
pdeplot(model,'XYData',sol(:,i,T(4)))
X2=[X 'dt*' num2str(T(4)) ')'];
title(X2)
colormap(jet); %caxis([amim amax]);
axis equal
end
end