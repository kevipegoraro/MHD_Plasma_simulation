function [X,Y]=Pcoordenada(R,Z) %pega matrizes em coordenadas cartesianas e tranforma em matrizes em pseutoridais
for j=1:length(Z(1,:))
    for i=1:length(Z(:,1))
        if round(Z(i,j),3)==0
            X(i,j)=(R(i,j)-0.37);
            Y(i,j)=0;
        else
            alfa = round(acot((R(i,j)-0.37)./Z(i,j)),4);
            X(i,j)=Z(i,j)./sin(alfa); 
            Y(i,j) =alfa;
        end
    end
    
end