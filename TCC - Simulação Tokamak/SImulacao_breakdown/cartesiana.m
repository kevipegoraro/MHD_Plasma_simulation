function Xi=cartesiana(X,Y,Z)
 if Z == 0 
     if sqrt(X^2+Y^2)==0.37
       r=0.01;
     end
 end
    %syms x
    Ng=3*65;
    x = linspace(0.0001,3*0.06,Ng);
    u = linspace(0,0,Ng);
    i=1; u(i) = round([0.37+x(i)*sqrt(1-(Z/x(i))^2)]*sqrt(1-(Y/(0.37+x(i)*sqrt(1-(Z/x(i))^2)))^2),3);
    mim = abs(X-u(i)); indice = 1;
    for i=1:Ng
         u(i) = round([0.37+x(i)*sqrt(1-(Z/x(i))^2)]*sqrt(1-(Y/(0.37+x(i)*sqrt(1-(Z/x(i))^2)))^2),3);
         aux = abs(X-u(i));
         if aux < mim
             mim = aux;
             indice = i;
         end
    %theta = arcsin(Z/x); 
    %fhi = arcsin(Y/(0.37+x*cos(theta)));
    %sqrt(1-(Z/x)^2)  (Y/(0.37+x*cos(arcsin(Z/x)))
    %eqn = X == [0.37+x*sqrt(1-(Z/x)^2)]*sqrt(1-(Y/(0.37+x*sqrt(1-(Z/x)^2)))^2) ;
    Xi(1) = x(indice);
    %Xi(1) = (vpasolve(eqn, x, [0.001 0.07]),4); % calculo r
    Xi(2) = real(round(asin(Z/Xi(1)),4));% calculo theta
    Xi(3) = real(round(asin(Y/(0.37+Xi(1)*cos(Xi(2)))),4)); % calculo fhi

end