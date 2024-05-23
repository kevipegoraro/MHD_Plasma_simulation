function build_greentable_nova

%% Defining grid size
pos_bv_in = [0.161 0.265];
pos_bv_ex = [0.53 0.34];
r=0.12/2; 
p_central = 0.37;
Ng = 65;
R1 = 0.06;
Z0 = 0;
dZ = 0.05;
NR = 4;

%% Clearing figure
figure(1)
clf
drawnow

%% Calculating grid points
out.r = linspace(p_central-3*R1,p_central+3*R1,Ng);
out.z = linspace(Z0-3*R1,Z0+3*R1,Ng);
[R,Z] = meshgrid(out.r,out.z);
%Total = zeros(size(R)).';
%% Loading machine coils
%W = load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/Adaptacao_nova/nova_vacuum_vessel');
%EF = load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/Adaptacao_nova/nova_coils');
%cnames = {'E1' 'E2' 'E3' 'E4'};
Ro = p_central; %0.1610;
Zo=0; %0.2650;
dR=0.01;
dZ=0.01;
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
            %disp(['Coil ' cnames{icoils} ': [R Z] =  [' num2str(ir) ' ' num2str(iz) ']'])
            G  = G  + 1/NR/NZ*green_em_nova(Rfilaments(ir),Zfilaments(iz),R,Z,'psi').';
            BR = BR + 1/NR/NZ*green_em_nova(Rfilaments(ir),Zfilaments(iz),R,Z,'br').';
            BZ = BZ + 1/NR/NZ*green_em_nova(Rfilaments(ir),Zfilaments(iz),R,Z,'bz').';
        end
    end
    out.G = G;  out.BR = BR;  out.BZ = BZ;
    % Plotting
    figure(1)
    clf
    %Total=G.'+Total;
    contour(out.r,out.z,G.',51,'b')
    hold on
    
    alfa = 0:0.1:(2*pi);
    xx = cos(alfa)*r;
    yy = sin(alfa)*r;
    plot(p_central+xx,yy,'r','linewidth',2)
    
    plot(pos_bv_in(1),pos_bv_in(2),'xr','linewidth',3)
    plot(pos_bv_in(1),-pos_bv_in(2),'xr','linewidth',3)
    plot(pos_bv_ex(1),pos_bv_ex(2),'xr','linewidth',3)
    plot(pos_bv_ex(1),-pos_bv_ex(2),'xr','linewidth',3)
    axis([p_central-4*r,p_central+4*r,-2*r,2*r]);
    
    xlabel('R (m)')
    ylabel('Z (m)')
    hold off
    axis equal
    grid on
    drawnow
    pause(2)
end
    %hold on
    %contour(out.r,out.z,Total,51,'b')
