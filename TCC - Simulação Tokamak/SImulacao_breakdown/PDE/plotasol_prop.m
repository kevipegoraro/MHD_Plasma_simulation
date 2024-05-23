function plotasol_prop(model,sol,c)
global G e
s=length(sol(1,1,:));
T=[3 round(s/4,0) round(3*s/4,0) s-1];
for i=c
    if i==1
        X=['n(:,'];
    elseif i==2
        X=['J_p_h_i(:,'];
    elseif i==3
        X=['Pe(:,'];
    elseif i==4
        X=['Pi(:,'];
    end
amax=max(max(sol(:,i,:)./sol(:,i,1)));
figure(i)
subplot(2,2,1)
pdeplot(model,'XYData',sol(:,i,T(1))./sol(:,i,1))
X2=[X 'dt*' num2str(T(1)) ')'];
title(X2)
colormap(jet); caxis([1 amax]);
axis equal 
subplot(2,2,2)
pdeplot(model,'XYData',sol(:,i,T(2))./sol(:,i,1))
X2=[X 'dt*' num2str(T(2)) ')'];
title(X2)
colormap(jet); caxis([1 amax]);
axis equal 
subplot(2,2,3)
pdeplot(model,'XYData',sol(:,i,T(3))./sol(:,i,1))
X2=[X 'dt*' num2str(T(3)) ')'];
title(X2)
colormap(jet); caxis([1 amax]);
axis equal 
subplot(2,2,4)
pdeplot(model,'XYData',sol(:,i,T(4))./sol(:,i,1))
X2=[X 'dt*' num2str(T(4)) ')'];
title(X2)
colormap(jet); caxis([1 amax]);
axis equal
end
end