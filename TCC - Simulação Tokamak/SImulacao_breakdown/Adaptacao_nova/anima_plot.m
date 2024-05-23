function anima_plot(corrente_plasma,corrente_bobinasparalelas,corrente_bobinas,ativa_quiver,m)
format short
close all
for i=1:1:m
        %escolhas das correntes %anima_plot(1,10,6,1,20)
    corrente_plasma = corrente_plasma*(1+0.4*(i/m)); %potencia do campo da corrente de plasma
    %for j=1:1:m
        j=0;
        corrente_bobinasparalelas = corrente_bobinasparalelas*(1+j/m); %potencia do campo da corrente de plasma
        %for k=1:1:m
        k=0;
             corrente_bobinas = corrente_bobinas*(1+0.1*(k/m));%potencia do campo da corrente de plasma   
    
             %chama a plotagem
             plot_circulo(corrente_plasma,corrente_bobinasparalelas,corrente_bobinas,ativa_quiver)
             pause(0.001)
        %end
    %end
end
end