n=10000;
A = zeros(1, n);
Asum = zeros(1, n);
Numberzeros = zeros(1, n);
Numberzeros(1)=1; Numberzeros(2)=2;
sum=0;
for i=3:1:n
    c=0;
    cont=0;
    d=0;
    while c==0 
        cont=cont+1;
        if A(round(i-1-cont,0)) == A(round(i-1,0))
            A(i)=cont;
            sum=sum+cont;
            c=1;
            d=1;
        end
        if round((i-cont),0)==2
            c=1;
        end
        if d==0
            Numberzeros(i)=Numberzeros(i-1)+1;
        else
            Numberzeros(i)=Numberzeros(i-1);
        end
    end
    Asum(i)=sum; %log(sum)

end
%%
x=1:1:n;
%plot(x,A,x,log(Asum));
figure(2)
prop0 = A./Numberzeros;
prop1 = A./Asum;
prop2 = Numberzeros./Asum;
prop3= zeros(1,n); for i=1:1:n-1 prop3(i)=A(i+1)/(1+A(i)); end
s=zeros(1,n); for i=2:1:n s(i)=(i*s(i-1)+prop3(i))/(i+1); end
plot(x,prop2,x,prop0,x,prop2,x,prop1,x,prop2,x,prop3)
%figure(3)
%plot(prop2,s)
%%