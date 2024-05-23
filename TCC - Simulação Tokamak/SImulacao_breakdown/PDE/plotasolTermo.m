function plotasolTermo(model,sol)
%figure(6); 
%pdemesh(model); 
%axis equal
%te=pe/n/e
e= 1.6e-19;
s=length(sol(1,1,:));
T=[3 round(s/4,0) round(3*s/4,0) s-1];
X=['T_e(:,'];
figure(9)
subplot(2,2,1)
aux=sol(:,3,T(1))./sol(:,1,T(1))/e;
pdeplot(model,'XYData',aux)
X2=[X 'dt*' num2str(T(1)) ')'];
title(X2)
colormap(jet);
axis equal 
subplot(2,2,2)
aux=sol(:,3,T(2))./sol(:,1,T(2))/e;
%ind = find(dn<n0);
%    dn(ind) = n0;
pdeplot(model,'XYData',aux)
X2=[X 'dt*' num2str(T(2)) ')'];
title(X2)
colormap(jet);
axis equal 
subplot(2,2,3)
aux=sol(:,3,T(3))./sol(:,1,T(3))/e;
pdeplot(model,'XYData',aux)
X2=[X 'dt*' num2str(T(3)) ')'];
title(X2)
colormap(jet);
axis equal 
subplot(2,2,4)
aux=sol(:,3,T(4))./sol(:,1,T(4))/e;
pdeplot(model,'XYData',aux)
X2=[X 'dt*' num2str(T(4)) ')'];
title(X2)
colormap(jet);
axis equal


X=['T_e_,_i(:,'];
figure(10)
subplot(2,2,1)
aux=(sol(:,4,T(1))+sol(:,3,T(1)))./sol(:,1,T(1))/e;
pdeplot(model,'XYData',aux)
X2=[X 'dt*' num2str(T(1)) ')'];
title(X2)
colormap(jet);
axis equal 
subplot(2,2,2)
aux=(sol(:,4,T(2))+sol(:,3,T(2)))./sol(:,1,T(2))/e;
%ind = find(dn<n0);
%    dn(ind) = n0;
pdeplot(model,'XYData',aux)
X2=[X 'dt*' num2str(T(2)) ')'];
title(X2)
colormap(jet);
axis equal 
subplot(2,2,3)
aux=(sol(:,4,T(3))+sol(:,3,T(3)))./sol(:,1,T(3))/e;
pdeplot(model,'XYData',aux)
X2=[X 'dt*' num2str(T(3)) ')'];
title(X2)
colormap(jet);
axis equal 
subplot(2,2,4)
aux=(sol(:,4,T(4))+sol(:,3,T(4)))./sol(:,1,T(4))/e;
pdeplot(model,'XYData',aux)
X2=[X 'dt*' num2str(T(4)) ')'];
title(X2)
colormap(jet);
axis equal
end