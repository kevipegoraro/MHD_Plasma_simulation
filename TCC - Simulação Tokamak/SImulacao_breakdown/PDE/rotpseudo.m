function [Dr, Dtheta, Dphi] = rotpseudo(fr,ftheta,fphi,dtheta,dr,raio,theta)
    R0=0.37;
    naux=size(fphi);
    nR=naux(1);
    ntheta=naux(2);
    Dr=zeros(nR,ntheta); Dtheta=Dr; Dphi=Dr; %raio(i)*R0+cos(theta(j))raio(i)^2
    for i=1:(nR-1)
    for j=1:(ntheta-1)
       deriR=(raio(i+1)*ftheta(i+1,j)-raio(i)*ftheta(i,j))/dr/raio(i); 
       deriTheta =(fr(i,j+1)-fr(i,j))/dtheta/raio(i); 
       Dphi(i,j)=deriR-deriTheta;
       Dr(i,j) =( (R0+cos(theta(j+1))*raio(i))*fphi(i,j+1)-(R0+cos(theta(j))*raio(i))*fphi(i,j))/dtheta/(raio(i)*R0+cos(theta(j))*raio(i)^2);
       Dtheta(i,j) =-( (R0+cos(theta(j))*raio(i+1))*fphi(i+1,j)-(R0+cos(theta(j))*raio(i))*fphi(i,j))/dr/(R0+cos(theta(j))*raio(i));
    end
    end
    for j=1:ntheta Dtheta(end,j) = 0;  end
    for i=1:nR     Dtheta(i,end)=Dtheta(i,1);  end
    
    for j=1:ntheta Dphi(end,j)=0;  end     
    for i=1:nR     Dphi(i,end)=Dphi(i,1);  end
    
    for j=1:ntheta Dr(end,j)=0;  end 
    for i=1:nR     Dphi(i,end)=Dphi(i,1);  end
end