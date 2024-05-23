function [Maa,Mvv,Mav,Mva,Raa,Rvv] = create_em_model

%% Loading TCABR coils
psiem = load('/home/canal/matlab/equil_reconst/greenfunc/green_table_65x65.mat');
psiemh = load('/home/canal/matlab/equil_reconst/greenfunc/green_table_129x129.mat');
Coil = load('/home/canal/matlab/equil_reconst/greenfunc/tcabr_coils');
C = Coil([1 8:23],:);
cname = {'OH' 'E1' 'E2' 'E3' 'E4' 'E5' 'E6' 'E7' 'E8' 'E9' 'E10' 'F1' 'F2' 'F3' 'F4' 'D1' 'D2'};

%% Active-active coil coupling
% Resistance Matrix
Raa = 1.68e-8*Coil(:,5)*2*pi.*Coil(:,1)./(Coil(:,3).*Coil(:,4)./Coil(:,5)); %Surface area divided by the number of turns
Raa = diag([sum(Raa(1:7)); Raa(8:23)]);

%Inductance matrix
L = self_inductance(C(:,1),C(:,3).*C(:,4),C(:,5));
L(1,1) = 0;
for ioh = 1:7
    L(1,1) = L(1,1) + max(interp2(psiemh.r,psiemh.z,Coil(ioh,5)*psiemh.OH.G.',linspace(Coil(ioh,1)*0.8,Coil(ioh,1)*1.2,1001),Coil(ioh,2)*ones(1,1001),'cubic'));
end
Maa = diag(L);
for icoil = 1:size(C,1)
    for jcoil = icoil+1:size(C,1)
        eval(['Maa(icoil,jcoil) = mutual_inductance(psiem.r,psiem.z,psiem.' cname{icoil} '.G,C(jcoil,1),C(jcoil,2),C(jcoil,3),C(jcoil,4),C(jcoil,5),10);'])
    end
end
Maa = triu(Maa) + triu(Maa,1)';

%% Vessel-vessel filament coupling
[W,dRv,dZv] = vessel_filaments;
dA = dRv.*dZv;

% Resistance matrix
Rvv = diag(7.4e-7*2*pi*W(:,1)./dA);

Mvv = zeros(size(W,1),size(W,1));
%Inductance matrix
for icoil = 1:size(W,1)
    Mvv(icoil,:) = green_em(W(icoil,1),W(icoil,2),W(:,1),W(:,2),'mutual');
    Mvv(icoil,icoil) = self_inductance(W(icoil,1),5e-5,1);
end

%% Active-vessel coupling
%Inductance matrix
Mav = zeros(size(C,1),size(W,1));
for icoil = 1:size(C,1)
    eval(['Mav(icoil,:) = interp2(psiem.r,psiem.z,psiem.' cname{icoil} '.G.'',W(:,1),W(:,2),''cubic'');']);
end
Mva = Mav.';

if 0

%% Main circuit matrices
M = [Maa Mav; Mva Mvv]; % Mutual inductances
R = blkdiag(Raa,Rvv); % Resistances

%% Plasma-plasma (needs the plasma current density)
%Inductance matrix
rp = linspace(0.42,0.825,65); drp = mean(diff(rp));
zp = linspace(-0.24,0.24,65); dzp = mean(diff(zp));
[zzp,rrp] = meshgrid(zp,rp);
rrp = rrp(:);
zzp = zzp(:);
Mpp = zeros(length(rrp));
for ip = 1:length(rrp)
    Mpp(ip,:) = green_em(rrp(ip),zzp(ip),rrp,zzp,'mutual');
    Mpp(ip,ip) = 0;
end

%% Active-plasma (needs the plasma current density)
%Inductance matrix
Map = zeros(size(C,1),length(rrp));
for icoil = 1:size(C,1)
    eval(['Map(icoil,:) = interp2(psiem.r,psiem.z,psiem.' cname{icoil} '.G.'',rrp,zzp,''cubic'');']);
end

%% Vessel-plasma
%Inductance matrix
Mvp = zeros(size(W,1),length(rrp));
for ivv = 1:size(W,1)
    Mvp(ivv,:) = green_em(W(ivv,1),W(ivv,2),rrp,zzp,'mutual');
end

end

end