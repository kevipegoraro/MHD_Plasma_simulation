function y1=massamolapendulo(x0,y0,h)
xf=0.5;
x = x0:h:15; 
y1 = zeros(4, length(x));
y1(4,1) = y0(4);
y1(3,1) = y0(3);
y1(2,1) = y0(2);
y1(1,1) = y0(1);
i = x0;
kevi = 1
for i=1:length(x)-1  %funcao der(valor do tempo, vetor com 4 entradas
        k_1 = der4(x(i),y1(:,i));
        k_2 = der4(x(i)+0.5*h,(y1(:,i)+0.5*h*k_1));
        k_3 = der4((x(i)+0.5*h),(y1(:,i)+0.5*h*k_2));
        k_4 = der4((x(i)+h),(y1(:,i)+k_3*h));
        y1(1,i)
        y1(:,i+1) = y1(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h; 
end

end