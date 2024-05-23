%Exemplo de gradiente em coordenadas pseudoroidais
%campo spim
%posições tokamak
r = 0.06;
p_central = 0.37;
Ng = 65;
%out.r = linspace(0.001,2*0.06,Ng); %out.theta = linspace(0,2*pi,Ng);
%%out.fhi = linspace(0,2*pi,Ng); %[R,Theta] = meshgrid(out.r,out.theta);
% criar a mash grid em coordnadas cartesianas
x = linspace(p_central-2*r,p_central+2*r,Ng);
y = linspace(-2*r,2*r,Ng);
[X,Z] = meshgrid(x,y);
hold off
%definindo a função pra trabalhar no campo de vetores
%F=X+Z;  %[DX,DY] = gradient((F),.24,.24);
DX=(X-Z); %./sqrt(X.^2+Z.^2); 
DZ=DX.^2+DZ.^2;
figure(1)
quiver(X,Z,DX,DZ,0.6) 
title('quiver em cartesianas');
axis([p_central-2*r p_central+2*r -2*r 2*r]); %espaço de posições em cartesianas
%%%%% fazendo a figura em pseutoridais
[b11, b12] = Pcoordenada(X,Z); %pega matrizes em coordenadas cartesianas e tranforma em matrizes em pseutoridais
[b21, b22] = Pcoordenada(DX,DZ); 
figure(2)
quiver(b11, b12,b21, b22,0.6) %vetores
axis([0 2*r -pi pi]); %espaço de posições em pseutoridais
title('espaço de posições em pseutoridais')
%%% fazendo a figura voltando pra cartesianas
[X1,Y1,Z1] = Pcoordenada_pseudo_to_cart(b11,b12,0); %pega matrezes de cordenadas pseutoridais e manda pra cartesianas
[X2,Y2,Z2] = Pcoordenada_pseudo_to_cart(b21,b22,0); %pega matrezes de cordenadas pseutoridais e manda pra cartesianas
figure(3)
quiver(X1,Z1,X2,Z2,0.6);
axis([p_central-2*r p_central+2*r -2*r 2*r]); %espaço de posições em cartesianas
title('espaço de posições em pseutoridais trnformadas em cartesianas')




%contour(X,Y,(F_plasma),50) %linhas de campo
%4 campos devido as bobinas 
%[r,theta,fhi] = meshgrid(out.r,out.theta,out.fhi);
% e=ones(Ng,Ng,Ng,3);
% e2=ones(Ng,Ng,Ng,3);
% for i = 1:Ng
%     for j = 1:Ng
%         for k = 1:Ng
%             a = pseudotoridal(out.r(i),out.theta(j),out.fhi(k));
%             e(i,j,k,1) = a(1); e(i,j,k,2) = a(2); e(i,j,k,3) = a(3);
%             n = sqrt(2*0.37*cos(out.theta(j))+.37^2+out.r(i)^2);
%             X = -e(i,j,k,2)/n; Y = e(i,j,k,1)/n; Z = e(i,j,k,1)/n;
%             b = cartesiana(X, Y, Z);
%             e2(i,j,k,1) = b(1); e2(i,j,k,2) = b(2); e2(i,j,k,3) = b(3);
%         end
%     end
% end

% e=ones(Ng,Ng,3);
% e2=ones(Ng,Ng,3);
% for i = 1:Ng
%     for j = 1:Ng
% %         for k = 1:Ng
%              a = pseudotoridal(out.r(i),out.theta(j),0);
%              e(i,j,1) = a(1); e(i,j,2) = a(2); e(i,j,3) = a(3);
%              n = sqrt(2*0.37*cos(out.theta(j))+.37^2+out.r(i)^2);
%              X = -a(2)/n; Y = a(1)/n; Z = a(3)/n;
%              b = cartesiana(X, Y, Z);
%              e2(i,j,1) = b(1); e2(i,j,2) = b(2); e2(i,j,3) = b(3);
% %         end
%      end
% end
%    a = pseudotoridal(rr,out.theta(i),0);

 
%%
% close all
% %clear all
% clc
% hold on
% for j=1:10
% for i=1:10
%       %   figure(1)
%          j=5;
%         plot([e(i,j,1) e(i,j,2)],[e2(i,j,1) e2(i,j,2)],'-'); axis([0 0.06 0 2*pi]);
%         %plot([k(j)-0.37 k2(j)],[-k2(j) k(j)-0.37],'-'); axis([-2*pi 2*pi -2*pi 2*pi]);
%        %      [cord x1   cord y1],[cord x2 cord y2]
%         %  figure(2)
%        %  plot(e(i,j,1), e2(i,j,1),'x'); axis([0 1 0 1]); 
% %plot([rr e2(i,1)],[out.theta(i) e2(i,2)],'-'); axis([0 1 0 2*pi]);
% %plot3([e(i,1) e2(i,1)],[e(i,2) e2(i,2)],[e(i,3) e2(i,3)],'-'); axis([0 1 0 2*pi 0 2*pi]);
% end
% end
%%