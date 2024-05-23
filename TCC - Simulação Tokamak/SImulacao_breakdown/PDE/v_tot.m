function f = v_tot(location,state)
e = 1.6e-19;
me = 9.11e-31;
mi = 1.67e-27;
%N = 1; nr = length(location.x); f = zeros(N,nr); 
% retornara: (v_ei + v_en + v_ion - v_loss)
f = 0;
