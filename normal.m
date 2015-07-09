function prob = normal(x,C);

invC = inv(C);
detC = det(C);

prob = exp(-0.5*x*invC*x');
prob = prob / sqrt(detC) / (2*pi)^(length(x)/2);

return
