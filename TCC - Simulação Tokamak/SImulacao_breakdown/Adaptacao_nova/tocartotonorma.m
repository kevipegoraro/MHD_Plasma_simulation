function a=tocartotonorma(v,l) %pega cordenadas pseutoridais e manda pra cartesianas
 X1=pseudotoridal2(v);
 X2=pseudotoridal2(l);
 a=norm(X2-X1); 
end