function F = integration(x,f,F0)
if length(x) == length(f)
    f(isnan(x) | isnan(f)) = 0;
    F = zeros(size(f)); F(1) = F0;
    for ii = 1:length(x) - 1
        F(ii+1) = F(ii) + ( f(ii) + f(ii+1) )*(x(ii+1) - x(ii))/2;
    end
else
    error('Vectors must have same length')
end
end