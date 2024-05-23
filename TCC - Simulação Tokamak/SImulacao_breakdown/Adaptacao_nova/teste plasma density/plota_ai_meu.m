function plota_ai_meu(n,R,Z)
        colormap(jet); 
        [q,h]=contourf(R,Z,n,40); 
        set(h,'linecolor','none');
        axis equal;
        xlabel('R (m)'); ylabel('Z (m)');
        colorbar; 
        %axis([-0.07 0.07 -0.07 0.07]);
 
end