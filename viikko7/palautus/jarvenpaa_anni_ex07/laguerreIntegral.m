function I = laguerreIntegral(N)

I = 0;
syms x;
laguerre = laguerreL(N, x);
r = double(solve(laguerre, x));

for i = 1:N
    I = I + r(i) / (( (N+1)^2) * (laguerreL(N+1, r(i)))^2) * cos(r(i))^2;
end

end
    