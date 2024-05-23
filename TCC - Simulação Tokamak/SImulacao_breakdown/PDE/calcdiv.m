function div=calcdiv(Jr,Jz,dr,dz)
div = zeros(size(Jr));
div1=div; div2=div;
for i=1:(length(Jr(1,:))-1)
div1(i,:) = (Jr(i+1,:)-Jr(i,:))/dr;
div2(:,i) = (Jz(:,i+1)-Jz(:,i))/dz;
end
div1(end,:) = (Jr(end,:)-Jr(end-1,:))/dr;
div2(:,end) = (Jz(:,end)-Jz(:,end-1))/dz;
div=div1+div2;
end