function [X,Y,Z]=Pcoordenada_pseudo_to_cart(R,theta,fhi) %pega matrezes de cordenadas pseutoridais e manda pra cartesianas
if length(fhi) > 1
for j=1:length(theta(1,:))
    for i=1:length(theta(:,1))
            o = pseudotoridal(R(i,j),theta(i,j),fhi(i,j)); %acot((R(i,j)-0.37)./Z(i,j));
            X(i,j)=o(1); 
            Y(i,j) =o(2);
            Z(i,j) = o(3);
        end
    end
else
    for j=1:length(theta(1,:))
    for i=1:length(theta(:,1))
            aux = (0.37+R(i,j)*cos(theta(i,j)));            
            X(i,j)= aux*round(cos(fhi),2);
            Y(i,j)= aux*round(sin(fhi),2);
            Z(i,j) = R(i,j)*sin(theta(i,j));
        end
    end
end
end