load('twophase_sampling.mat'); 
m = size(uu, 1); dt = tt(2) - tt(1); 
nodes = [0, logspace(-1, 1, 9)/4, inf];
g = {}; 
for nodecount = 1:length(nodes)-1
    g{nodecount} = util_bandpass(nodes(nodecount), nodes(nodecount+1), 2);
end
z = zeros(m, length(g)+1, length(tt));
for i = 1:m
    z(i, 1, :) = yy(i, :);
    for gcount = 1:length(g)
        z(i, gcount+1, :) = lsim(g{gcount}, uu(i, :), tt);
    end
end
gam = zeros(m, length(g)+1, length(g)+1);
for i = 1:m
    for k1 = 1:size(gam, 2)
        for k2 = 1:size(gam, 3)
            temp1 = reshape(z(i, k1, 1:end), 1, length(tt));
            temp2 = reshape(z(i, k2, 1:end), length(tt), 1);
            gam(i, k1, k2) = dt*(temp1*temp2);
        end
    end
end

omegas = logspace(-2, 2, 101)/4; 
nu4 = 0.05; leg4 = ['Butterworth ($\nu=', num2str(nu4, '%.2f'), '$)'];
[M4, rho4, xi4] = util_ocsvm(gam, nu4);
psi_found4 = util_frequency_evaluate(M4, g, omegas); 
nu5 = 0.20; leg5 = ['Butterworth ($\nu=', num2str(nu5, '%.2f'), '$)'];
[M5, rho5, xi5] = util_ocsvm(gam, nu5);
psi_found5 = util_frequency_evaluate(M5, g, omegas); 
nu6 = 0.40; leg6 = ['Butterworth ($\nu=', num2str(nu6, '%.2f'), '$)'];
[M6, rho6, xi6] = util_ocsvm(gam, nu6);
psi_found6 = util_frequency_evaluate(M6, g, omegas); 

fprintf('Butterworth. Softness = %.2f \t Average violation = %.4e \t Margin = %.4e \t\n', nu4, mean(xi4), rho4);
fprintf('Butterworth. Softness = %.2f \t Average violation = %.4e \t Margin = %.4e \t\n', nu5, mean(xi5), rho5);
fprintf('Butterworth. Softness = %.2f \t Average violation = %.4e \t Margin = %.4e \t\n', nu6, mean(xi6), rho6);

figure; legends = {leg1, leg2, leg3, leg4, leg5, leg6}; 
semilogx(omegas*4, psi_found1, 'Color', [209 73 5]/255);
hold on, semilogx(omegas*4, psi_found2, 'Color', [153 0 0]/255);
hold on, semilogx(omegas*4, psi_found3, '--m');
hold on, semilogx(omegas*4, psi_found4, 'Color', [0 132 115]/255);
hold on, semilogx(omegas*4, psi_found5, 'Color', [65 86 161]/255);
hold on, semilogx(omegas*4, psi_found6, '--c');
legend(legends, 'interpreter', 'latex'); 
xlabel('$\omega\tau_0$', 'Interpreter', 'latex'),
ylabel('$\ell(j\omega)$', 'Interpreter', 'latex'); 