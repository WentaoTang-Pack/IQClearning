%% Step response examination
stepsim = @(t, x) twophase_model(x, [1; 0; 0], [0; 0]); 
t_span = 0:0.1:45;
[~, sol] = ode45(stepsim, t_span, zeros(6, 1)); 
y_act = sol(:, 2); 
sys_nom = tf(0.28, [4, 1], 'InputDelay', 12);
y_nom = lsim(sys_nom, 0*t_span + 1, t_span);
figure; 
plot(t_span, y_act, 'k', t_span, y_nom, '--k');
xlabel('$t$ (min)', 'Interpreter', 'latex'), ylabel('$y(t)$', 'Interpreter', 'latex'), 
ylim([-0.05, 0.35]), legend({'Actual', 'Nominal'}); 