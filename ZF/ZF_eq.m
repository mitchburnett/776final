% ZF Equalizer
clearvars;

%% Load in an observation model -- Options ['Ungerboeck', 'Forney', 'PSMFOM']
model = 'Ungerboeck';

load(['../Models/', model, '.mat']);
L1 = eval([model, '.L1']);
L2 = eval([model, '.L2']);
hn = eval([model, '.filter']);
if strcmp(model, 'Ungerboeck')
    noise_filter = eval([model, '.noise_filter']);
end
eval(['clear ' model]);
%%

% Simulation data and parameters
Nsymbols = 10e4;
A = 1;

LUT = A*[-1-1j;...
         -1+1j;...
         1-1j;...
         1+1j];
  
idx = randi([1,4],Nsymbols,1);
I = LUT(idx);

SNR = 15; %dB
EbNo = 10.^(SNR./10); % snr natural number
var_n = A^2/(2*EbNo) * sum(abs(hn).^2) ;
std_n = sqrt(var_n);

s = conv(hn, I);
s = s(1:end-L2+1); % I still have a question as to why we dont cut off the front. Since a real data signal would not have come yet.

n = std_n*(randn(length(s),1) + 1j*randn(length(s),1));

% create correlated noise for Ungerboeck
if strcmp(model, 'Ungerboeck')
    n = filter(noise_filter, 1, n);
end

vn = s + n;

%%
% ZF equalizer
P = 40;  % Equalizer length.
n0 = 15; % Needs to be >= L1
[w_zf, err] = spike(hn, n0, P);

% Apply filter
yn = conv(w_zf, vn);
yn = yn(n0+1:end-P+n0+L1);
% So I spent some time trying to get correct where the data symbols and
% the estimate line up. I am confued why we can't cut off the anti causal
% part in the first convolution.

hn_eq = conv(w_zf,hn);

%%
nfft = 1024;
H = fftshift(fft(hn, nfft));
H_ZF = fftshift(fft(w_zf, nfft));
fspan = linspace(-.5, .5, nfft);

%% Plots
figure(1); clf;
plot(vn,'.'); grid on; hold on;
plot(yn, '.');
axis('square');

% Look at the ZF EQ filter.
figure(2);
subplot(311)
stem(L1:L2, abs(hn)); grid on;
title('|h(nT)|');
subplot(312);
stem(0:P-1,abs(w_zf)); grid on;
title('FIR Least Squares approx. h(nT)');

subplot(313);
stem(0:length(hn_eq)-1, abs(hn_eq)); grid on;

figure(3); clf;
plot(fspan, 20*log10(abs(H))); grid on; hold on;
plot(fspan, 20*log10(abs(H_ZF)));
legend('H', 'ZF');
