function [X,Y]=Pcoordenada_contour(Z) %pega matrizes em coordenadas cartesianas e tranforma em matrizes em pseutoridais
for j=1:length(Z(:,1)) %z( ,j) represneta valor de Z
    for i=1:length(Z(1,:)) %z(i,) representa valor de R
        if round(Z(i,j),3)==0
            X(i,j)=(Z(i,j)-0.37);
            Y(i,j)=0;
        else
            alfa = round(acot((R(i,j)-0.37)./Z(i,j)),4);
            X(i,j)=Z(i,j)./sin(alfa); 
            Y(i,j) =alfa;
        end
    end
    
end