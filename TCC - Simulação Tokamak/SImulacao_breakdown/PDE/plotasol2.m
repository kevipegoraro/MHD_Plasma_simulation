function plotasol2(model,sol)
%figure(6); 
%pdemesh(model); 
%axis equal

%TEMPOS PARA PLOTAR
T=[1 4 10 20];
%keyboard
    %DENSIDADE DE PARTICULAS PLOTS
subplot(4,4,3)
pdeplot(model,'XYData',sol(:,1,T(1)))
X=['n(:,']; X2=[X 'dt*' num2str(T(1)) ')']; title(X2); colormap(jet); axis([-1 1 -1 1]); 

subplot(4,4,4)
pdeplot(model,'XYData',sol(:,1,T(2)))
X=['n(:,']; X2=[X 'dt*' num2str(T(2)) ')']; title(X2); colormap(jet); axis([-1 1 -1 1]); 

subplot(4,4,5)
pdeplot(model,'XYData',sol(:,1,T(3)))
X=['n(:,']; X2=[X 'dt*' num2str(T(3)) ')']; title(X2); colormap(jet); axis([-1 1 -1 1]);  

subplot(4,4,6)
pdeplot(model,'XYData',sol(:,1,T(4)))
X=['n(:,']; X2=[X 'dt*' num2str(T(4)) ')']; title(X2); colormap(jet); axis([-1 1 -1 1]);


%DENSIDADE DE CORRENTE PLOTS

subplot(4,4,7)
pdeplot(model,'XYData',sol(:,2,T(1)))
X=['J_p_h_i(:,']; X2=[X 'dt*' num2str(T(1)) ')']; title(X2); colormap(jet); axis([-1 1 -1 1]);


subplot(4,4,8)
pdeplot(model,'XYData',sol(:,2,T(2)))
X=['J_p_h_i(:,']; X2=[X 'dt*' num2str(T(2)) ')']; title(X2); colormap(jet); axis([-1 1 -1 1]);


subplot(4,4,9)
pdeplot(model,'XYData',sol(:,2,T(3)))
X=['J_p_h_i(:,']; X2=[X 'dt*' num2str(T(3)) ')']; title(X2); colormap(jet); axis([-1 1 -1 1]);  

subplot(4,4,10)
pdeplot(model,'XYData',sol(:,2,T(4)))
X=['J_p_h_i(:,']; X2=[X 'dt*' num2str(T(4)) ')']; title(X2); colormap(jet); axis([-1 1 -1 1]);


%PLOTAGEM PRESSAO ELETRONICA
subplot(4,4,11)
pdeplot(model,'XYData',sol(:,3,T(1)))
X=['Pe(:,']; X2=[X 'dt*' num2str(T(1)) ')']; title(X2);  colormap(jet); axis([-1 1 -1 1]);

subplot(4,4,12)
pdeplot(model,'XYData',sol(:,3,T(2)))
X=['Pe(:,']; X2=[X 'dt*' num2str(T(2)) ')']; title(X2); colormap(jet); axis([-1 1 -1 1]);

subplot(4,4,13)
pdeplot(model,'XYData',sol(:,3,T(3)))
X=['Pe(:,']; X2=[X 'dt*' num2str(T(3)) ')']; title(X2); colormap(jet); axis([-1 1 -1 1]);

subplot(4,4,14)
pdeplot(model,'XYData',sol(:,3,T(4)))
X=['Pe(:,']; X2=[X 'dt*' num2str(T(4)) ')']; title(X2); colormap(jet); axis([-1 1 -1 1]);

%PLOAGEM PRESSAO IONICA
subplot(4,4,15)
pdeplot(model,'XYData',sol(:,4,T(1)))
X=['Pi(:,']; X2=[X 'dt*' num2str(T(1)) ')']; title(X2); colormap(jet); axis([-1 1 -1 1]);

subplot(4,4,16)
pdeplot(model,'XYData',sol(:,4,T(4)))
X=['Pi(:,']; X2=[X 'dt*' num2str(T(4)) ')']; title(X2); colormap(jet); axis([-1 1 -1 1]);
end