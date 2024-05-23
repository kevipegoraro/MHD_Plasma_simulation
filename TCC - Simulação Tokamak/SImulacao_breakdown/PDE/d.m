function c = d(location, state)
%keyboard
c = ones(4,length(location.x));
c(1,:)=location.x;
c(2,:)=0.0001;
end