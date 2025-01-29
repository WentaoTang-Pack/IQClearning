% load('delay_sampling.mat'); 
m = size(uu, 1); dt = tt(2) - tt(1); ndelay = round(delay/dt); 
g0 = tf([0, 1], [1/pi^2, sqrt(2)/pi, 1]); 
g1 = tf([1/pi^2, 0, 0], [1/pi^2, sqrt(2)/pi, 1]); g2 = g1 * g0;
g = [g0, g1, g2]; 
z = zeros(m, length(g)+1, length(tt));
for i = 1:m
    z(i, 1, :) = yy(i, :);
    for gcount = 1:length(g)
        z(i, gcount+1, :) = lsim(g(gcount), uu(i, :), tt);
    end
end
gam = zeros(m, length(g)+1, length(g)+1);
for i = 1:m
    for k1 = 1:size(gam, 2)
        for k2 = 1:size(gam, 3)
            temp1 = reshape(z(i, k1, ndelay+1:end), 1, length(tt)-ndelay);
            temp2 = reshape(z(i, k2, ndelay+1:end), length(tt)-ndelay, 1);
            gam(i, k1, k2) = dt*(temp1*temp2);
        end
    end
end

omegas = logspace(-2, 2, 101); 
[M2, rho2, xi2] = util_ocsvm(gam, 0.05, m, 3);
psi_found2 = util_frequency_evaluate(M1, [g0, g1, g2], omegas);