%função condições iniciais
function uinit = inifuncadapt_adapt(location)
global Neq Te0 n0 e 
%keyboard
r=1-(location.x.^2+location.y.^2).^2;
uinit = zeros(Neq,1); %cria com todos zeroslength(location.x)
uinit(1,:) = r*n0;
uinit(2,:)=-r*1e-6;
uinit(3,:)=r*e*Te0*n0;
uinit(4,:)=r*e*Te0*n0;
