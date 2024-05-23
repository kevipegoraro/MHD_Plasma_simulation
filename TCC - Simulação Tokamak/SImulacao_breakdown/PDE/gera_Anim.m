function gera_Anim(model,sol,animoquem)
global G e
if G.plotton
    for oll=2:(G.cplott-1)
        eval(['auxx = G.plott' num2str(oll) '(1);' ])
        if isnan(auxx) ==0
            figure(11);
            X=['(vin-vlos)(:,' num2str(oll) ')'];
            %figure(11+oll)
            eval(['pdeplot(model,''XYData'',G.plott' num2str(oll) ');' ])
            colormap(jet); title(X);
           % keyboard
            axis equal
            figure(12);
            X=['Campoplasma(:,' num2str(oll) ')'];
            %figure(11+oll)
            eval(['pdeplot(model,''XYData'',G.plottcampoplasma' num2str(oll) ');' ])
            colormap(jet); title(X);
            axis equal
            pause(0.3)
        end
    end
end

if animoquem(1) %se animaquem indice 1 for 1, terei a animacao de uma propriedade da solucao
    figure(12);
    s=length(sol(3,2,:));
    %p=1;
    %if s>50 p=3; end
    %if s>100 p=8; else 
    p=round(s*0.15/10,0) 
    %end 
    T=1:p:(s-1); %T da os tempos que iremos ver na animacao 
    i=animoquem(2); %animaquem indice 2 dis qual a propriedade da solucao que sera animada
        if i==1
            X=['n(:,'];
        elseif i==2
            X=['J_p_h_i(:,'];
        elseif i==3
            X=['Pe(:,'];
        elseif i==4
            X=['Pi(:,'];
        end
        amax=max(max(-sol(:,i,:)));
        amim=max(max(sol(:,i,:)));
    for oll=T
        pdeplot(model,'XYData',-sol(:,i,oll))
        X2=[X 'dt*' num2str(oll) ')'];
        title(X2)
        colormap(jet);
        axis equal
        %caxis([amim amax])
        pause(0.15)
    end
end