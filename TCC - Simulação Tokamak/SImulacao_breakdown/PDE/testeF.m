systemsize = 10;
location.x = 0;
location.y = 0;  
location.z = 0;  
location.subdomain=1;  
state.u = zeros(systemsize, 1);  
state.ux = zeros(systemsize, 1);  
state.uy = zeros(systemsize, 1);  
state.uz = zeros(systemsize, 1);  
state.time = 0;
% Now test the function call:  
mhandle=@s03;  
result = mhandle(location, state)