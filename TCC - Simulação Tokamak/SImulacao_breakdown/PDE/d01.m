function c = d01(location, state)
c = ones(4,length(location.x));
c(1,:)=location.x;
c(2,:)=0.001;
