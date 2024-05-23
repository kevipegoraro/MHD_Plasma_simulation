%função condições iniciais
function uinit = inifunc_meio(location,state)
%primeira edp U= zeros(Neq,1); U(1)=n0; U(5)=e*Te0*n0; U(6)=U(5);
global Neq Te0 n0 e p0 
%keyboard 
aj = 1e-6*(1-sqrt(locaxion.x(1)^2+locaxion.y(1)^2))
uinit = zeros(Neq,length(location.x)); %cria com todos zeros
uinit(1,:) = n0;
uinit(2,:)= aj;
uinit(3,:)=aj;
uinit(4,:)=1e-6;
uinit(5,:)=e*Te0*n0;
uinit(6,:)=e*Te0*n0;
uinit(7,:)= 7.89e11*p0; %v_en 
%uinit(8,:)= 0; %v_ion 
uinit(9,:) = 0; %out1.v_ei; %v_ei
%uinit(10,:)= 0; %v_loss
