load('delay_sampling.mat'); 
m = size(uu, 1); dt = tt(2) - tt(1); ndelay = round(delay/dt); 
g0 = tf([0, 1], [1, 1]); g1 = tf([1, 0], [1, 1]); g2 = g1 * g0;
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
[M1, rho1, xi1] = util_ocsvm(gam, 0.05, m, 3);
psi_found1 = util_frequency_evaluate(M1, [g0, g1, g2], omegas);


% [M2, rho2, xi2] = ocsvm(gam, 0.2, m, 3);
% psi_found2 = util_frequency_evaluate(M2, [g0, g1, g2], omegas);
% [M3, rho3, xi3] = ocsvm(gam, 0.8, m, 3);
% psi_found3 = util_frequency_evaluate(M3, [g0, g1, g2], omegas);

psi_theor = 4*sin(omegas*delay/2).^2 .* (omegas<pi/delay) + 4*(omegas>=pi/delay); 
psi_megre = (4*omegas.^4 + 50*omegas.^2) ./ (omegas.^4 + 6.5*omegas.^2 + 50);
figure; semilogx(omegas, psi_theor, 'k', omegas, psi_megre, '--k');
hold on, semilogx(omegas, psi_found1, 'Color', [0.8500 0.3250 0.0980]);
hold on, semilogx(omegas, psi_found2, 'm');
% hold on, semilogx(omegas, psi_found3, 'r');
legends = {'Actual', 'Megretski-Rantzer', 'Learned', 'Learned (Butterworth)'};
legend(legends, 'interpreter', 'latex'); 
xlabel('$\omega$', 'Interpreter', 'latex'), ylabel('$\ell(j\omega)$', 'Interpreter', 'latex'); 