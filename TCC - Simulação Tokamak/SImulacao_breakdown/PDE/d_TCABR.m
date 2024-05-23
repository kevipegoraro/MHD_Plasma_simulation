function c = d_TCABR(location, ~)
c = ones(4,length(location.x));
c(1,:)=location.x;
c(2,:)=0;
