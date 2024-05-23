function a=pseudotoridal2(v) %pega cordenadas pseutoridais e manda pra cartesianas
r=v(1);
theta=v(2);
fhi=v(3);
aux = round((0.37+r*cos(theta)),9);
a(1) = round(aux*cos(fhi),9); % calculo do X
a(2) = round(aux*sin(fhi),9); % calcul odo Y
a(3) = round(r*sin(theta),9);
end