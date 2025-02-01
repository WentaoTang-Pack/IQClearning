A = 1; T = 2000; dt = 0.2; tt = 0:dt:T; 
sys_nom = tf(0.28, [4, 1], 'InputDelay', 12);

figure; 
colors = lines(5); 

omegas = 7.5*[1e-3, 1e-2, 1e-1, 1e0]; 
for j = 1:length(omegas)
    omega = omegas(j); 
    sim = @(t, x) twophase_model(x, [A*cos(omega*t); 0; 0], [0; 0]); 
    tt = 0:dt:60;
    [~, sol] = ode45(sim, tt, zeros(6, 1)); 
    response_actual = sol(:, 2)';
    response_nominal = lsim(sys_nom, A*cos(omega*tt), tt);
    hold on, subplot(2, 2, j), plot(tt, response_actual, 'Color', colors(j, :)); 
    hold on, subplot(2, 2, j), plot(tt, response_nominal, 'Color', colors(j, :), 'LineStyle', '--'); 
    legends = cell(1, 2); 
    legends{1} = ['Actual ($\omega\tau_0 =' , num2str(omega*4), '$)']; 
    legends{2} = ['Nominal ($\omega\tau_0 =' , num2str(omega*4), '$)']; 
    hold on, subplot(2, 2, j), legend(legends, 'interpreter', 'latex'), ylim([-0.3, 0.305]), 
    xlabel('$t$ (min)', 'Interpreter', 'latex'), ylabel('$y(t)$', 'Interpreter', 'latex');
end

