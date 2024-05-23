
Ng1=40; Ng2=Ng1*2;
Ntime = 50; %1000 ou 10000
dt = 1e-5; %tempo correto é 1e-4
R0 = 0.37;
B0 = 0.7;
Vloop = 10;  %corrente gerando o campo E externo
p0 = 0.05; %0.01 % 1e-4 milibar = 0.01 Pa
ng = p0/1.38e-23/300; %particulas do gás
n0 = 1e2; %densidade de particulas inicial do plasma
Te0 = 0.026; % teperatura inicial do gas
gas = 'H2'; %tipo de gas
a0=0.07;
dR = a0/(Ng1+1); dZ = 2*pi/Ng2; %Definindo dr e dz
r_ini=dR;
raio = linspace(r_ini,a0,Ng1);
theta = linspace(0,2*pi,Ng2);
[R,Z] = meshgrid(raio,theta);
Rcart = R.*cos(Z);
Zcart = R.*sin(Z);

n=n0*ones(Ng1,Ng2,Ntime+1); %definindo a condição inicial da densidade de particulas
x0=a0/2;
y0=0;
sgma=-1;
s0=1;
D=0.001;
po=1;
%s=zeros(Ng1,Ng2);
s=s0*exp(((+a0/2.5+Rcart').^2+(Zcart').^2)/sgma);
%plota_ai_meu(s',Rcart,Zcart)
% keyboard
for it = 1:Ntime
    laplaciano = D*Polar_laplaciano(n(:,:,it),dR,dZ,raio,n0); 
    figure(2); plota_ai_meu(laplaciano,Rcart,Zcart)
    for aux1 = 0:(Ng1-1)
        n((end-aux1),:,it+1) = n((end-aux1),:,it) + dt*(s((end-aux1),:)+laplaciano((end-aux1),:));
    end
     n(:,1,it+1) = n(:,end,it+1); n(:,2,it+1) = n(:,end-1,it+1); %valores em theta 0 são igual a theta 2pi
                                   % %n(1,:,it+1) = n(1,1,it);
    if po==1  
        figure(1);
        plota_ai_meu(n(:,:,it+1),Rcart,Zcart); ti1=['n(:,:,' num2str(it) ')'];  title(ti1);
        po=0;
       keyboard
    end
    po=po+1;
end
    