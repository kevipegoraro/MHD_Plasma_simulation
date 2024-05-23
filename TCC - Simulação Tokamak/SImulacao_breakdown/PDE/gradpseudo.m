function [Dr, Dtheta] = gradpseudo(f,dtheta,dr,raio)
    naux=size(f);
    nR=naux(1);
    ntheta=naux(2);
    Dr=zeros(nR,ntheta); Dtheta=Dr;
    for i=1:(nR-1)
    for j=1:(ntheta-1)
       Dr(i,j)=(f(i+1,j)-f(i,j))/dr; %/(raio(i));
       Dtheta(i,j) = (f(i,j+1)-f(i,j))/(dtheta*raio(i));
    end
    end
    for j=1:ntheta    Dtheta(i,end) = (f(i,end)-f(i,end-1))/(dtheta*raio(i));   end
    for i=1:(nR)  Dr(end,j)=(f(end,j)-f(end-1,j))/dr; end
end