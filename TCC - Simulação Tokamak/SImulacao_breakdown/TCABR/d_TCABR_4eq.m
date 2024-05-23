function c = d_TCABR_4eq(location, ~)
c = zeros(4,length(location.x));
c(1,:)=location.x;
c(2,:)=9.11e-31;
c(3,:)=3/2;
c(4,:)=3/2;
