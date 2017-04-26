  
Tfsamp = 10; % samples/symbol
figC =1;
SPAN = 12;
SL = -SPAN/2;
SR = SPAN/2;
% pulse shape
g = rcosdesign(.5, 12, 10).'; % 50% excess bandwidth, spanning 12 symbols, 10 samples/symbol.
gn = g(1:Tfsamp:end);
g_span = linspace(SL, SR, length(g));
gn_span = linspace(SL, SR, length(gn)); 

% plot of the pulse shape and downsampled pulse shape gn = g(nT).
figure(figC); figC=figC+1; clf;
plot(g_span,g); grid on; hold on;
stem(gn_span,gn); hold off;
title('pulse shape - g(t)');
legend('g(t)', 'gn = g(nT)');

% channel model
c = zeros(Tfsamp+1,1);
c(1) = 1;
c(Tfsamp+1) = -.9;

% plot the channel model.
figure(figC); figC=figC+1;
stem((0:length(c)-1)/Tfsamp, c); grid on;
xlabel('t/T (symbols)');
title('c(t)');

% h(t) = g(t)*c(t);
h = conv(g,c);
h_span = linspace(-6,7,length(h));

% plot h(t)
figure(figC); figC=figC+1;
plot(h_span, abs(h)); grid on;
title('|h(t)|');

% x(t) = h(t)*h(-t)
x = conv(h,flip(h));
x_span = linspace(-13,13,length(x));

% plot x(t)
figure(figC); figC=figC+1;
plot(x_span,abs(x)); grid on;
title('|x(t)|');

% compute frequency domain for each signal.
nfft = 1024;
G = fftshift(fft(g,nfft));
C = fftshift(fft(c,nfft));
H = fftshift(fft(h,nfft));
X = fftshift(fft(x,nfft));

% plot G(f) and C(f)
figure(figC); figC=figC+1;
freq_span = linspace(-5,5,nfft);
plot(freq_span, 20*log10(abs(G))); grid on; hold on;
plot(freq_span, 20*log10(abs(C))); grid on; hold off;
ylim([-40, 20]);
legend('G(f)', 'C(f)');
ylabel('Magnitude (dB)');

% plot G(f) and H(f)
figure(figC); figC=figC+1;
plot(freq_span, 20*log10(abs(G))); grid on; hold on;
plot(freq_span, 20*log10(abs(H))); grid on; hold off;
ylim([-40, 20]);
xlim([-5, 5]);
legend('G(f)', 'H(f)');
ylabel('Magnitude (dB)');

% plot G(f) and X(f)
figure(figC); figC=figC+1;
plot(freq_span, 20*log10(abs(G))); grid on; hold on;
plot(freq_span, 20*log10(abs(X))); grid on; hold off;
ylim([-20, 40]);
xlim([-5, 5]);
legend('G(f)', 'X(f)');
ylabel('Magnitude (dB)');

% downsampled xn = x(nT)
xn = x(1:Tfsamp:end);
figure(figC); clf; figC=figC+1;
plot(x_span,abs(x)); grid on; hold on;
xn_span = linspace(-13,13, length(xn));
plot(xn_span, abs(xn), 'o'); hold off;
title('x(t) and xn=x(nT)');
legend('|x(t)|', 'xn');

% shorten xn to 99.99% energy
xn_short = xn((abs(xn) > .1));
figure(figC); clf; figC=figC+1;
stem([-1, 0, 1], abs(xn_short)); grid on;
title('Shorten xn');
legend('|xn|');

% compute DTFT of xn
XN = fftshift(fft(xn_short,nfft));
xn_span = linspace(-.5, .5, nfft);
figure(figC); clf; figC=figC+1;
plot(xn_span, 20*log10(abs(XN))); grid on;
ylim([-50, 50]);
title('DTFT of xn');
xlabel('frequency (cycles/symbol)');
ylabel('Magnitude (dB)');

% perform spectral factorization
xn_roots = roots(xn_short);
figure(figC); clf; figC=figC+1;
zplane(xn_roots); grid on;
title('Spectral Factorization of X(z)');

f_n = [1, -xn_roots(2)];

% normalize f_n to unit energy.
sigma_sq = sum(abs(f_n).^2);
sigma = sqrt(sigma_sq);

f_n = f_n/sigma;

% plot |f_n|
figure(figC); clf; figC=figC+1;
stem([0, 1], abs(f_n)); grid on;
title('|f_n|');

% compute DTFT of F(z)
F_z = fftshift(fft(f_n, nfft));
F_z_span = linspace(-0.5, 0.5, nfft);

% plot F(z)
figure(figC); clf; figC=figC+1;
plot(F_z_span, 20*log10(abs(F_z))); grid on;
ylim([-40, 20]);
title('DTFT of F(z)');
ylabel('Magnitude (dB)');


