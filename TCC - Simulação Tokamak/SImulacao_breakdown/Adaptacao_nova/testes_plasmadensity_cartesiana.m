clc;
clf;
Ntime = 50; %1000 ou 10000
dt = 1e-2; %tempo correto é 1e-4
n0 = 1e2; %densidade de particulas inicial do plasma
%%%%%%%%%%
Ng1=65; Ng2=65;
a0=0.06;
R0=0.37;
ri = linspace(-a0,a0,Ng1);
zi = linspace(-a0,a0,Ng2);
[R,Z] = meshgrid(ri,zi);
dR = 2*a0/Ng1; %mean(diff(R)); 
dZ = 2*a0/Ng2; %mean(diff(Z)); 


%%%%
n= n0*ones(Ng1,Ng2,Ntime+1); %definindo a condição inicial da densidade de particulas
x0=-a0/3;
y0=0.03;
sgma=1;
s0=0.1;
s=s0*exp((2-(R+x0)^2+(Z+y0)^2)/sgma);
D=0.0001;
 %keyboard
%%%%%%%%%%%
po=0;
for it = 1:Ntime
    laplaciano = D*del2(n(:,:,it),dR,dZ);
    figure(2); 
    plota_ai_meu(laplaciano,R,Z)
    n(:,:,it+1) = n(:,:,it) + dt*(s+laplaciano);
    n(1,:,it+1) = n0; 
    n(end,:,it+1)=n0;
    n(:,end,it+1) = n0; 
    n(:,1,it+1) = n0;
    if po==2
    figure(1);
    colormap(jet); [q,h]=contourf(R,Z,n(:,:,it),Ng1); set(h,'linecolor','none');
    axis equal;
    xlabel('R (m)'); ylabel('Z (m)');
    colorbar; axis([-0.07 0.07 -0.07 0.07]);
    ti1=['n(:,:,' num2str(it) ')'];
    title(ti1);
    po=0;
    keyboard
    end
    po=po+1;
end
 