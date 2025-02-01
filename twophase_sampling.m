%% Sampling
A = 0.25; T = 45; dt = 0.1; tt = 0:dt:T;  
m = 500; % Number of trajectories
uu = zeros(m, length(tt)); yy = zeros(m, length(tt));
freqs = sort(10.^(-2 + 4*rand(1, m))/4); 
sys_nom = tf(0.28, [4, 1], 'InputDelay', 12);
for i = 1:m 
    uu(i, :) = A*sin(freqs(i)*tt); 
    % actual
    % forward Euler
    % current_x = zeros(6, 1);
    % response_actual = zeros(length(tt), 1);
    % for tcount = 1:length(tt)-1
    %     v = twophase_model(current_x, [uu(i, tcount);0;0], [0;0]);
    %     current_x = current_x + v*dt;
    %     response_actual(tcount+1) = current_x(2);
    % end
    % Runge-Kutta
    sim = @(t, x) twophase_model(x, [A*sin(freqs(i)*t); 0; 0], [0; 0]); 
    [~, sol] = ode45(sim, tt, zeros(6, 1)); 
    response_actual = sol(:, 2); 
    % nominal
    response_nominal = lsim(sys_nom, uu(i, :), tt);
    yy(i, :) = response_actual - response_nominal;
end
figure; colors = spring(25);
for j = 1:25
    hold on, subplot(1, 2, 1), plot(tt, uu(20*j, :), 'Color', colors(j, :)), 
    xlabel('$t$', 'Interpreter', 'latex'), ylabel('$u(t)$', 'Interpreter', 'latex');
    hold on, subplot(1, 2, 2), plot(tt, yy(20*j, :), 'Color', colors(j, :)), 
    xlabel('$t$', 'Interpreter', 'latex'), ylabel('$y(t)$', 'Interpreter', 'latex');
end
save('twophase_sampling.mat', 'uu', 'yy', 'tt', 'A', 'T', 'freqs'); 
savefig('twophase_sampling.fig');