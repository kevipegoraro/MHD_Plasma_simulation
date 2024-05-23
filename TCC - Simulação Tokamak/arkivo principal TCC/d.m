function c = d(location, state)
c = ones(4,length(location.x));
c(1,:)=location.y;
c(2,:)=0.00001;
