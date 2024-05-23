function build_greentable_nova_unidadetrocada
clear all 
close all
clc
%% Defining grid size
pos_bv_in = [0.161 0.265];
pos_bv_ex = [0.53 0.34];
r11=0.06; 
p_central = 0.37;

isave = 1;
Ng = 65;
R1 = 0.06;
NR = 2;

%% Calculating grid points
%gerar mash grid retangulares
out.r = linspace(p_central-2*R1,p_central+2*R1,Ng);
out.z = linspace(-3*R1,3*R1,Ng);
[R,Z] = meshgrid(out.r,out.z);
%gerar mash grid polares
out.r1 = linspace(0.001,2*r11,Ng); % r represena a distancia até o ponto central
out.z1 = linspace(0,pi*(2-2/Ng),Ng); %z representa o angulo alfha n minah cordenada "polar"
[R1,Z1] = meshgrid(out.r,out.z);

%% Loading machine coils
%W = load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/Adaptacao_nova/nova_vacuum_vessel');
EF = load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/Adaptacao_nova/nova_coils');
cnames = {'E1' 'E2' 'E3' 'E4' 'E50' 'E51' 'E52' 'E53' 'E6' 'E7' 'E8' 'E9'};



%% Calculating Green's functions between external coils and grid points
for icoils = 1:length(cnames)

    Ro = EF(icoils,1);
    Zo = EF(icoils,2);
    dR = EF(icoils,3);
    dZ = EF(icoils,4);
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
            %disp(['Coil ' cnames{icoils} ': [R Z] =  [' num2str(ir) ' ' num2str(iz) ']'])
            G  = G  + 1/NR/NZ*green_em_nova_unidadestroca(Rfilaments(ir),Zfilaments(iz),R,Z,'psi').';
            BR = BR + 1/NR/NZ*green_em_nova_unidadestroca(Rfilaments(ir),Zfilaments(iz),R,Z,'br').';
            BZ = BZ + 1/NR/NZ*green_em_nova_unidadestroca(Rfilaments(ir),Zfilaments(iz),R,Z,'bz').';
        end
    end
    eval(['out.' cnames{icoils} '.G = G;'])
    eval(['out.' cnames{icoils} '.BR = BR;'])
    eval(['out.' cnames{icoils} '.BZ = BZ;'])
    
end

    X=Pcoordenada(EF(icoils,1),EF(icoils,2))
    Ro1 = X(1);
    Zo1 = X(2);


if isave
    G = out;
    fname = ['/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/green_table/green_table_nova_coredanadaspolares_' num2str(Ng) 'x' num2str(Ng) '.mat'];
    save(fname,'-struct','G')
    disp(['Green''s table saved at ' fname])
end


end