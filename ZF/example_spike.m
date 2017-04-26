% example spike

%g(n) = -0.2delta(n) + delta(n-2)
g = [-0.2, 0, 1];

n = 13;
n0 = 8;

% compute the least squares inverse approximate.
[h, err] = spike(g, n0, n);

% compute the equalized signal
d = conv(g,h);

g_plot = zeros(1,n);
g_plot(1) = g(1);
g_plot(3) = g(3);

figure(99);
subplot(311)
stem(0:n-1, g_plot); grid on;
title('g(n)');
ylim([-.5, 1.25]);
subplot(312);
stem(0:n-1,h); grid on;
title('FIR Least Squares approx. h(n)');
ylim([-.5, 1.25]);
subplot(313);
stem(0:n-1,d(1:n)); grid on;
title('Equalized signal for delay n_0=8')
ylim([-.5, 1.25]);

