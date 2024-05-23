function aux = Polar_laplaciano(n,dR,dZ,raio)
%\nabla^2 f = \frac{\partial^2 f}{\partial \rho^2} + \frac{1}{\rho} 
%\frac{\partial f}{\partial \rho} + \frac{\partial^2 f}{\partial z^2}
naux=size(n);
nR=naux(1);
nZ=naux(2);
a=zeros(nR+2,nZ); %colocando dois pontos virtuais a mais de raio
soma=0;
for j=1:nZ
    soma=soma+n(1,j);
end
a(1,:)=soma/(nZ); %fazendo o ponto central virtual valer a média dos pontos imediatamnte apos ele.
a(end,:)=0;%n0;
b=a;
c=a;
n2=n;
n=a;
n(2:(end-1),:)=n2(:,:);
%size(a)
%size(n)
for i=2:(nR)
    for j=2:(nZ-1)
        a(i,j) = (n(i+1,j)-n(i-1,j))/dR/(raio(i));
        %derivada primeira ordem com relação a r
        b(i,j) = (n(i+1,j)-2*n(i,j)+n(i-1,j))/(dR^2);
        %derivada segunda com relação a r, dR
        c(i,j) = (n(i,j+1)-2*n(i,j)+n(i,j-1))/(dZ^2)/(raio(i)^2);
        %derivada segunda ordem com relaçã oa oangulo
    end
end
aux=a(2:(end-1),:)+b(2:(end-1),:)+c(2:(end-1),:);
aux=aux';
end