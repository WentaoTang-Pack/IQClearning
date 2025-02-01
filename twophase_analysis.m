load('twophase_sampling.mat'); 
m = size(uu, 1); dt = tt(2) - tt(1); 
nodes = [0, logspace(-1, 1, 9)/4, inf];
g = {}; 
for nodecount = 1:length(nodes)-1
    g{nodecount} = util_bandpass(nodes(nodecount), nodes(nodecount+1));
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
nu1 = 0.05; leg1 = ['Simple ($\nu=', num2str(nu1, '%.2f'), '$)'];
[M1, rho1, xi1] = util_ocsvm(gam, nu1);
psi_found1 = util_frequency_evaluate(M1, g, omegas); 
nu2 = 0.1; leg2 = ['Simple ($\nu=', num2str(nu2, '%.2f'), '$)'];
[M2, rho2, xi2] = util_ocsvm(gam, nu2);
psi_found2 = util_frequency_evaluate(M2, g, omegas); 
nu3 = 0.2; leg3 = ['Simple ($\nu=', num2str(nu3, '%.2f'), '$)'];
[M3, rho3, xi3] = util_ocsvm(gam, nu3);
psi_found3 = util_frequency_evaluate(M3, g, omegas); 

fprintf('Softness = %.2f \t Average violation = %.4e \t Margin = %.4e \t\n', nu1, mean(xi1), rho1);
fprintf('Softness = %.2f \t Average violation = %.4e \t Margin = %.4e \t\n', nu2, mean(xi2), rho2);
fprintf('Softness = %.2f \t Average violation = %.4e \t Margin = %.4e \t\n', nu3, mean(xi3), rho3);

figure; legends = {leg1, leg2, leg3}; 
semilogx(omegas*4, psi_found1, 'Color', [209 73 5]/255);
hold on, semilogx(omegas*4, psi_found2, 'Color', [153 0 0]/255);
hold on, semilogx(omegas*4, psi_found3, 'Color', [250 200 0]/255);
legend(legends, 'interpreter', 'latex'); 
xlabel('$\omega\tau_0$', 'Interpreter', 'latex'),
ylabel('$\ell(j\omega)$', 'Interpreter', 'latex'); 

