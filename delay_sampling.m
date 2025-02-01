A = 1; K = 5; T = 5;
tt = 0:0.01:T; delay = 0.5;
ndelay = round(delay/(tt(2)-tt(1))); 
m = 500; % Number of trajectories
uu = zeros(m, length(tt)); yy = zeros(m, length(tt));
freqs = zeros(1, m);
for i = 1:m 
    % Frequency sampling
    freqs(i) = 10^(-2 + 4*rand); 
    uu(i, :) = A*sin(freqs(i)*tt); 
    % Binary sequence 
    % uu(i, :) = A*sign(2*rand(1, length(tt)) - 1); 
    % Fourier coefficient sampling
    % Fourier_random_sampling(A, K, tt);
    % Output
    yy(i, :) = [zeros(1, ndelay), uu(i, 1:end-ndelay)] - uu(i, :);
end
figure;
subplot(1, 2, 1), plot(tt, uu(1:20:m, :)), 
xlabel('$t$', 'Interpreter', 'latex'), ylabel('$u(t)$', 'Interpreter', 'latex');
subplot(1, 2, 2), plot(tt, yy(1:20:m, :)), 
xlabel('$t$', 'Interpreter', 'latex'), ylabel('$y(t)$', 'Interpreter', 'latex');

save('delay_sampling.mat', 'uu', 'yy', 'tt', 'A', 'K', 'T', 'delay', 'freqs'); 
savefig('delay_sampling.fig');

function uu = Fourier_random_sampling(A, K, tt)
a0 = 0; % A*(2*rand-1); 
b = zeros(1, 5); b(1) = (2*rand - 1)*A;
for k = 2:K
    b(k) = (2*rand - 1)*abs(b(1))/k;
end
uu = zeros(1, length(tt)) + a0/sqrt(2);
for k = 1:K
    uu = uu + b(k)*sin(pi/2*k*tt/tt(end));
end
end
