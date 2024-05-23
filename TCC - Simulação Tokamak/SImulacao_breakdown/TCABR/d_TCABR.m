function c = d_TCABR(location, ~)
global G
c = zeros(6,length(location.x));
c(1,:)=location.x;
% ----teste----
c(2,:)=G.uxcoeficientes;
c(4,:)=G.uxcoeficientes;
% ----teste----
c(3,:)=9.11e-31;
c(5,:)=3/2;
c(6,:)=3/2;
