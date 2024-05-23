function div = divpseudo(fr,ftheta,dtheta,dr,raio,theta)
    R0=0.37;
    naux=size(f);
    nR=naux(1);
    ntheta=naux(2);
    Dr=zeros(nR,ntheta); Dtheta=Dr; %raio(i)*R0+cos(theta(j))raio(i)^2
    for i=1:(nR-1)
    for j=1:(ntheta-1)
       Dr(i,j)=( (raio(i-1)*R0+cos(theta(j))*raio(i+1)^2)*fr(i+1,j)-(raio(i)*R0+cos(theta(j))*raio(i)^2)*fr(i,j))/dr/(raio(i)*R0+cos(theta(j))*raio(i)^2); 
       Dtheta(i,j) =( (R0+cos(theta(j+1))*raio(i))*ftheta(i,j+1)-(R0+cos(theta(j))*raio(i))*ftheta(i,j))/dtheta/(raio(i)*R0+cos(theta(j))*raio(i)^2);
    end
    end
    for j=1:ntheta Dtheta(i,end) = Dtheta(i,1);  end
    for i=1:nR Dr(end,j)=0;end
    div=Dr+Dtheta; 
end