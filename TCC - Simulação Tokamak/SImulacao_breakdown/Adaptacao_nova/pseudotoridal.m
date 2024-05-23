function a=pseudotoridal(r,theta,fhi) %pega cordenadas pseutoridais e manda pra cartesianas
aux = (0.37+r*cos(theta));
a(1) = aux*cos(fhi);
a(2) = aux*sin(fhi);
a(3) = r*sin(theta);
end