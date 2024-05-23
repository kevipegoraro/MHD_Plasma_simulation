function plot_nova(fig)

%% Loading machine coils
%W = load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/Adaptacao_nova/nova_vacuum_vessel');
EF = load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/Adaptacao_nova/nova_coils');
% FL = load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/Adaptacao_nova/nova_flux_loops');
% MP = load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/Adaptacao_nova/nova_magnetic_probes');
% SC = load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/Adaptacao_nova/nova_saddle_coils');
pos_bv_in = [0.161 0.265];
pos_bv_ex = [0.53 0.34];
r=0.12/2; 
p_central = 0.37;
cnames = {'E1' 'E2' 'E3' 'E4' 'Ep0' 'Ep1' 'Ep2' 'Ep3' 'Ep4' 'E6' 'E7' 'E8'  'E9' 'E10'};
tam =15;%length(cnames)+1
%% Plotting
figure(fig)
hold on

% E-Coils
for iecoils = 1:tam
    plot([EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2 EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2 EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2],[EF(iecoils,2)-EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)-EF(iecoils,4)/2 EF(iecoils,2)-EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)-EF(iecoils,4)/2],'k')
end
% 
% % F-Coils
% for ifcoils = 8:size(EF,1)
%     plot([EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2 EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2 EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2],[EF(ifcoils,2)-EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)-EF(ifcoils,4)/2 EF(ifcoils,2)-EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)-EF(ifcoils,4)/2],'r')
% end

    r=0.06;
    prp=0.07/r;
    alfa = 0:0.1:(2*pi+0.1);
    xx = cos(alfa)*r;
    yy = sin(alfa)*r;
    plot(p_central+xx,yy,'r','linewidth',2)
    plot(p_central+xx*prp,yy*prp,'r','linewidth',2)
% Magnetic probes
if 0
    for improbe = 1:18
        improb = [1:9 17:25];
        imp = improb(improbe);
        plot([MP(imp,1)+2.5e-3 MP(imp,1)-2.5e-3 MP(imp,1)-2.5e-3 MP(imp,1)+2.5e-3 MP(imp,1)+2.5e-3],[MP(imp,2)-7.5e-3 MP(imp,2)-7.5e-3 MP(imp,2)+7.5e-3 MP(imp,2)+7.5e-3 MP(imp,2)-7.5e-3],'k')
    end
    for improbe = 1:14
        improb = [10:16 26:32];
        imp = improb(improbe);
        plot([MP(imp,1)+7.5e-3 MP(imp,1)-7.5e-3 MP(imp,1)-7.5e-3 MP(imp,1)+7.5e-3 MP(imp,1)+7.5e-3],[MP(imp,2)-2.5e-3 MP(imp,2)-2.5e-3 MP(imp,2)+2.5e-3 MP(imp,2)+2.5e-3 MP(imp,2)-2.5e-3],'k')
    end
end
% 
% %Flux loops
% if 0
%     plot(FL(:,1),FL(:,2),'kx','markersize',8,'linewidth',2)
% end
% 
% %Saddle coils
% if 0
%     %for isc = 1:size(SC,1)-2
%     %    plot([SC(isc,1)-0.035/2 SC(isc,1)+0.035/2],[SC(isc,2) SC(isc,2)],'-','color',[57 83 164]/256,'linewidth',3)
%     %end
%     plot([SC(19,1) SC(19,1)],[SC(19,2)-0.1 SC(19,2)+0.1],'o-','color',[57 83 164]/256,'linewidth',3)
% end
axis([p_central-3.5*r,p_central+3.5*r,-7*r,7*r]);
axis equal 
xlabel('R (m)')
ylabel('Z (m)')
grid on

end