function build_greentable_teste

%% Defining grid size
isave = 0;
Ng = 65;
R1 = 0.4;
R2 = 0.845;
Z0 = 0;
dZ = 0.26;
NR = 2;

%% Clearing figure
figure(1)
clf
drawnow

%% Calculating grid points
out.r = linspace(R1,R2,Ng);
out.z = linspace(Z0-dZ,Z0+dZ,Ng);
[R,Z] = meshgrid(out.r,out.z);

%% Loading machine coils
W = load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/tcabr_vacuum_vessel');
EF = load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/teste_coils');
cnames = {'E1' 'E2' 'E3' 'E4' 'E5' 'E6' 'E7' 'F1' 'F2' 'F3' 'F4' 'F5' 'F6' 'F7' 'F8' 'F9' 'F10' 'F11' 'F12' 'F13'};

%% Calculating Green's functions between external coils and grid points
for icoils = 1:4
    Ro = EF(icoils,1);
    Zo = EF(icoils,2);
    dR = EF(icoils,3);
    dZ = EF(icoils,4);
    NZ = round(NR*dZ/dR);
    Rin = Ro - dR/2 + dR/NR/2;
    Rout = Ro + dR/2 - dR/NR/2;
    Zin = Zo - dZ/2 + dZ/NZ/2;
    Zout = Zo + dZ/2 - dZ/NZ/2;
    Rfilaments = linspace(Rin,Rout,NR);
    Zfilaments = linspace(Zin,Zout,NZ);

    G = zeros(size(R)).'; BR = G; BZ = G;
    for ir = 1:NR
        for iz = 1:NZ
            disp(['Coil ' cnames{icoils} ': [R Z] =  [' num2str(ir) ' ' num2str(iz) ']'])
            G  = G  + 1/NR/NZ*green_em(Rfilaments(ir),Zfilaments(iz),R,Z,'psi').';
            BR = BR + 1/NR/NZ*green_em(Rfilaments(ir),Zfilaments(iz),R,Z,'br').';
            BZ = BZ + 1/NR/NZ*green_em(Rfilaments(ir),Zfilaments(iz),R,Z,'bz').';
        end
    end
    eval(['out.' cnames{icoils} '.G = G;'])
    eval(['out.' cnames{icoils} '.BR = BR;'])
    eval(['out.' cnames{icoils} '.BZ = BZ;'])

    % Plotting
    figure(1)
    clf
    contour(out.r,out.z,G.',51,'b')
    hold on
    for ir = 1:NR
        for iz = 1:NZ
            plot(Rfilaments(ir),Zfilaments(iz),'ok')
        end
    end
    for iecoils = 1:7
        plot([EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2 EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2 EF(iecoils,1)-EF(iecoils,3)/2 EF(iecoils,1)+EF(iecoils,3)/2],[EF(iecoils,2)-EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)-EF(iecoils,4)/2 EF(iecoils,2)-EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)+EF(iecoils,4)/2 EF(iecoils,2)-EF(iecoils,4)/2],'k')
    end
    for ifcoils = 8:size(EF,1)
        plot([EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2 EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2 EF(ifcoils,1)-EF(ifcoils,3)/2 EF(ifcoils,1)+EF(ifcoils,3)/2],[EF(ifcoils,2)-EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)-EF(ifcoils,4)/2 EF(ifcoils,2)-EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)+EF(ifcoils,4)/2 EF(ifcoils,2)-EF(ifcoils,4)/2],'r')
    end
    plot(W(:,1),W(:,2),'k')
    xlabel('R (m)')
    ylabel('Z (m)')
    hold off
    axis equal
    grid on
    axis([0.25 1.3 -0.699999 0.700001])
    drawnow
    pause(1)
end

if isave
    G = out;
    fname = ['/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/green_table/green_table_' num2str(Ng) 'x' num2str(Ng) '.mat'];
    save(fname,'-struct','G')
    disp(['Green''s table saved at ' fname])
end

end