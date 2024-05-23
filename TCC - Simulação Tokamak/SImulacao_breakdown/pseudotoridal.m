function Xi=pseudotoridal(r,theta,fhi) %pega cordenadas pseutoridais e manda pra cartesianas
aux = (0.37+r*cos(theta));
Xi(1) = aux*cos(fhi);
Xi(2) = aux*sin(fhi);
Xi(3) = r*sin(theta);
end