mus = 0.01:0.01:.99;
x1 = zeros(length(mus), 1);
x2 = zeros(length(mus), 1);
x3 = zeros(length(mus), 1);
syms x;
for i = 1:length(mus)
    eqn1 = x - (1-mus(i))/(x+mus(i))^2 - mus(i)/(x-(1-mus(i)))^2 == 0;
    eqn2 = x - (1-mus(i))/(x+mus(i))^2 + mus(i)/(x-(1-mus(i)))^2 == 0;
    eqn3 = x + (1-mus(i))/(x+mus(i))^2 - mus(i)/(x-(1-mus(i)))^2 == 0;
    x1(i) = vpasolve(eqn1, x, [1-mus(i), Inf]);
    x2(i) = vpasolve(eqn3, x, [-mus(i), 1-mus(i)]);
    x3(i) = vpasolve(eqn3, x, [-Inf, -mus(i)]);
end