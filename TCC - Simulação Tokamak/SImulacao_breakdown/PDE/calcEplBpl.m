function Mag=calcEplBpl(Jr,Jz,Jphi,n,dr,dz)
%%% A_pl -> A_r + A_z + A_phi 
%%% Modulo de A_pl = sqrt(A_r.^2 + A_z.^2 + A_phi.^2)
 % calculo de um divergente div=calcdiv(Jr,Jz,dr,dz)
 mi0=4*pi*1e-7; epsolon0=1/mi0/300000000^2;
A_r = zeros(size(Jr)); tam_r=length(A_r(:,1)); tam_z=length(A_r(1,:));
A_z=A_r; Aphi=A_r;
%soma das derivadas parciais do A_pl é igual a -mi0*J
for i=1:(tam_r-1)
A_r(i+1,:)= A_r(i,:)+dr*(-mi0*Jr(i,:));
end
A_r(end,:)= A_r(end-1,:)+dr*(-mi0*Jr(end,:));

for i=1:(tam_z-1)
A_z(:,i+1)= A_z(:,i)+dz*(-mi0*Jz(:,i));
end
A_z(:,end)= A_z(:,end-1)+dz*(-mi0*Jz(:,end));


A_phi(:,:,it+1)= A_phi(:,:,it)+dt*(-mi0*Jphi(:,:,it));


end