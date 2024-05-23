function L = self_inductance(R0,A,N)

L = 4e-7*pi*N.^2.*R0.*(log(8*R0./sqrt(A/pi)) - 7/4);

end