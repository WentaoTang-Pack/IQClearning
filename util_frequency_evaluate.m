function psi = util_frequency_evaluate(M, g, omegas)
j = sqrt(-1);
psi = zeros(1, length(omegas));
for k = 1:length(omegas)
    omega = omegas(k); 
    gval = zeros(length(g), 1); 
    for gcount = 1:length(g)
        gval(gcount) = polyval(g(gcount).Numerator{1}, j*omega)/polyval(g(gcount).Denominator{1}, j*omega); 
    end
    psi(k) = real(gval'*M*gval);
end