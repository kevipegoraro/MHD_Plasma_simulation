%somaNvezes(a, b, N) retorta a+b repetido N vezes
function aux = somaNvezes(a,b,N)
    aux=a+b;
    %x = 1:N
    for i=1:N,
        aux=aux+(a+b);
        %x[i-1]=aux
    end
    %plot(1:N,x)
    
    