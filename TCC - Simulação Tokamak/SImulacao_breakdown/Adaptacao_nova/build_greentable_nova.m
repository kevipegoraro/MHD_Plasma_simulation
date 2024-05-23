function build_greentable_nova
clear all 
close all
clc
%% Defining grid size
pos_bv_in = [0.161 0.265];
pos_bv_ex = [0.53 0.34];
r=0.12/2; 
p_central = 0.37;

isave = 1;
Ng = 65;
R1 = r;
NR = 2;

%% Calculating grid points
out.r = linspace(p_central-R1,p_central+R1,Ng);
out.z = linspace(-R1,R1,Ng);
%out.r = linspace(p_central-2*R1,p_central+2*R1,Ng);
%out.z = linspace(-3*R1,3*R1,Ng);
[R,Z] = meshgrid(out.r,out.z);

%% Loading machine coils
%W = load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/Adaptacao_nova/nova_vacuum_vessel');
EF = load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/Adaptacao_nova/nova_coils');
cnames = {'E1' 'E2' 'E3' 'E4' 'Ep0' 'Ep1' 'Ep2' 'Ep3' 'Ep4' 'E5' 'E6' 'E7' 'E8'  'E9' 'E10'};
t=length(cnames);
%% Calculating Green's functions between external coils and grid points
for icoils = 1:length(cnames)
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

    G = zeros(size(R)); BR = G; BZ = G;
    for ir = 1:NR
        for iz = 1:NZ
            %disp(['Coil ' cnames{icoils} ': [R Z] =  [' num2str(ir) ' ' num2str(iz) ']'])
            G  = G  + 1/NR/NZ*green_em_nova(Rfilaments(ir),Zfilaments(iz),R,Z,'psi');
            BR = BR + 1/NR/NZ*green_em_nova(Rfilaments(ir),Zfilaments(iz),R,Z,'br');
            BZ = BZ + 1/NR/NZ*green_em_nova(Rfilaments(ir),Zfilaments(iz),R,Z,'bz');
        end
    end
    eval(['out.' cnames{icoils} '.G = G'])
    eval(['out.' cnames{icoils} '.BR = BR;'])
    eval(['out.' cnames{icoils} '.BZ = BZ;'])
    
    % Plotting
    %figure(1)
    %clf
    %Total=G.'+Total;
    %contour(out.r,out.z,Total,51,'b')
    %hold on
    
%     alfa = 0:0.1:(2*pi);
%     xx = cos(alfa)*r;
%     yy = sin(alfa)*r;
%     plot(p_central+xx,yy,'r','linewidth',2)
%     
%     plot(pos_bv_in(1),pos_bv_in(2),'xr','linewidth',3)
%     plot(pos_bv_in(1),-pos_bv_in(2),'xr','linewidth',3)
%     plot(pos_bv_ex(1),pos_bv_ex(2),'xr','linewidth',3)
%     plot(pos_bv_ex(1),-pos_bv_ex(2),'xr','linewidth',3)
%     axis([p_central-4*r,p_central+4*r,-2*r,2*r]);
    
%     xlabel('R (m)')
%     ylabel('Z (m)')
%     hold off
%     axis equal
%     grid on
%     drawnow
%     pause(2)
end
    %keyboard
    %hold on
    %contour(out.r,out.z,Total,51,'b')
    
if isave
    G = out;
    fname = ['/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulacao_breakdown/green_table/green_table_nova_' num2str(Ng) 'x' num2str(Ng) '_2.mat'];
    save(fname,'-struct','G')
    disp(['Green''s table saved at ' fname])
end


%p = (out.E1.G+out.E2.G)*500+(out.E3.G+out.E4.G)*4000-out.E5.G*30000;
%figure(1); clf; contour(out.r,out.z,p.',51); axis equal
%p = (out.E1.G+out.E2.G)*500-(out.E3.G+out.E4.G)*500+out.E5.G*30000; colormap(jet); figure(1); clf; contour(out.r,out.z,p.',51); axis equal; plot_nova(1); colorbar;
%g = log(8*0.37/0.05)+0.05+1.2/2-1.5
%bz = 10^(-7)*30000*g/0.37;
end