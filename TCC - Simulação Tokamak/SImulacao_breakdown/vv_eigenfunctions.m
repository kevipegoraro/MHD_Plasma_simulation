function vv_eigenfunctions(Neig)

%% Getting Mvv Matrix and vacuum vessel wall
VV = load('/home/kevi/Área de Trabalho/TCC - Simulação Tokamak/SImulação_breakdown/tcabr_vacuum_vessel');
[~,Mvv,~,~,~,Rvv] = create_em_model;
[W,dRv,dZv] = vessel_filaments;

%% Vessel eigenmodes
A = -inv(Mvv)*Rvv;
[V,D] = eig(A);
[eigvalues,ind] = sort(diag(D),'descend');
V = V(:,ind(Neig));

%% Plotting
col = b2r_colormap(min(V),max(V));
dIv = (max(V) - min(V))/(size(col,1)-1);

figure(1)
clf
plot_tcabr(1)
for iv = 1:length(V)
    icol = 1 + floor((V(iv) - min(V))/dIv);
    fill([W(iv,1)-dRv(iv)/2 W(iv,1)+dRv(iv)/2 W(iv,1)+dRv(iv)/2 W(iv,1)-dRv(iv)/2 W(iv,1)-dRv(iv)/2],[W(iv,2)-dZv(iv)/2 W(iv,2)-dZv(iv)/2 W(iv,2)+dZv(iv)/2 W(iv,2)+dZv(iv)/2 W(iv,2)-dZv(iv)/2],col(icol,:),'edgecolor','none')
end
plot(VV(:,1),VV(:,2),'k')
colormap(b2r_colormap(min(V),max(V)))
colorbar
hold off
axis([0.2 1.3 -0.7 0.7])
title(['Eigenmode Decay Time: ' sprintf('%0.2f',-1/eigvalues(Neig)*1000) ' ms'])

figure(2)
clf
semilogy(eigvalues,'xk','linewidth',2)
hold on
semilogy(Neig,eigvalues(Neig),'xr','linewidth',3,'markersize',12)
hold off
title('Vacuum Vessel Eigenvalues ( s^{-1} )')
xlabel('Eigenvalue #')
axis([0 200 -10^5 -10^2])

figure(3)
clf
plot(1:size(V,1),V,'k','linewidth',2)
title(['Vacuum Vessel Eigenmode #' num2str(Neig)])
xlabel('Vacuum Vessel Filament #')

end