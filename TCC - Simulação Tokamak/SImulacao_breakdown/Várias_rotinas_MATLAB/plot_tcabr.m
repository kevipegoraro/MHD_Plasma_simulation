function plot_tcabr(fig)

%% Loading machine coils
W = load('/home/canal/matlab/machines/tcabr_vacuum_vessel');
EF = load('/home/canal/matlab/equil_reconst/greenfunc/tcabr_coils');
FL = load('/home/canal/matlab/equil_reconst/tcabr_flux_loops');
MP = load('/home/canal/matlab/equil_reconst/tcabr_magnetic_probes');
SC = load('/home/canal/matlab/equil_reconst/tcabr_saddle_coils');

%% Plotting
figure(fig)
hold on

% E-Coils
for iecoils = 1:7
    plot([EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2 EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2 EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2],[EF(iecoils,2)-EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)-EF(iecoils,4)/2 EF(iecoils,2)-EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)-EF(iecoils,4)/2],'k')
end

% F-Coils
for ifcoils = 8:size(EF,1)
    plot([EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2 EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2 EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2],[EF(ifcoils,2)-EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)-EF(ifcoils,4)/2 EF(ifcoils,2)-EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)-EF(ifcoils,4)/2],'r')
end

% Vacuum vessel
fill([W(1:5,1); W(7:11,1); W(1,1)],[W(1:5,2); W(7:11,2); W(1,2)],[1 1 1]*0.4,'edgecolor','none')
fill([W(7:11,1); W(13:17,1); W(7,1)],[W(7:11,2); W(13:17,2); W(7,2)],[1 1 1]*0.7,'edgecolor','none')
plot(W(:,1),W(:,2),'k')

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

%Flux loops
if 0
    plot(FL(:,1),FL(:,2),'kx','markersize',8,'linewidth',2)
end

%Saddle coils
if 0
    %for isc = 1:size(SC,1)-2
    %    plot([SC(isc,1)-0.035/2 SC(isc,1)+0.035/2],[SC(isc,2) SC(isc,2)],'-','color',[57 83 164]/256,'linewidth',3)
    %end
    plot([SC(19,1) SC(19,1)],[SC(19,2)-0.1 SC(19,2)+0.1],'o-','color',[57 83 164]/256,'linewidth',3)
end

xlabel('R (m)')
ylabel('Z (m)')
grid on
axis equal
axis([0.2 1.3 -0.7 0.7])

end