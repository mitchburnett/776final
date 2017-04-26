% MMSE Equalizer
%clearvars;
% load in an observation model.
load('../Models/Forney.mat');
fn = Forney.filter;
L = length(fn)-1;

% set up symbols
Nsymbols = 10e3;
A = 1;

LUT = A*[-1-1j;...
         -1+1j;...
         1-1j;...
         1+1j];
  
idx = randi([1,4],Nsymbols,1);
I = LUT(idx); % data

SNR = 15; %dB
EbNo = 10.^(SNR./10); % snr natural number

sigma_n = sqrt(A^2/(2*EbNo));

n = sigma_n*(randn(Nsymbols,1) + 1j*randn(Nsymbols,1));

r = I + n;

IF = conv(fn, I);

% I was doing IF = IF(length(f_n):end), however there were some floating symbols in the middle. Instead I am cutting the last L-1
%IF = IF(1:end-length(f_n)+1);
IF = IF(1:end-(length(fn)-1));
vn = IF + n;

figure(1); clf;
plot(I, '.'); grid on; hold on;
plot(r, '.');
%plot(IF, '.');
plot(vn,'.');

%% MMSE equalizer.

% create weiner hopf equations. This section I drew the approach from my
% 777 homework on matlab assignment c72.

% Gamma matrix
ADD = 1; % added this variable to experiment with the filter length that is most acceptable.
FN = convm(fn,2*L+ADD);
Gamma = zeros(2*L+1, 2*L+ADD);
for l = 1:2*L+ADD
    for i = 1:2*L+ADD
        Gamma(l,i) = FN(:,i).'*conj(FN(:,l)); % FN(:,i).'*conj(FN(:,l)) this seems to be the correct conjugate form.
    end
end

% Noise matrix
N = convm(n, 2*L+ADD);
Rn = zeros(2*L+1, 2*L+ADD);
for l = 1:2*L+ADD
    for i = 1:2*L+ADD
        Rn(l,i) = N(:,i)'*N(:,l)/(Nsymbols);
    end
end

Gamma = Gamma + Rn;

zeta = conj([flip(fn); zeros(L+ADD-1,1)]);

copt = Gamma\zeta;
fn_mmse = copt;

%%
nfft = 1024;
F = fftshift(fft(fn, nfft));
F_MMSE = fftshift(fft(fn_mmse, nfft));
fspan = linspace(-.5, .5, nfft);

figure(3); clf;
plot(fspan, 20*log10(abs(F))); grid on; hold on;
plot(fspan, 20*log10(abs(F_MMSE)));
legend('F', 'MMSE');

%%
yn = conv(fn_mmse, vn);
f_eq = conv(fn_mmse, fn);
%%
figure(1); clf;
%plot(yn(16:end-14), '.');
plot(I, '.'); grid on; hold on;
plot(r, '.');
plot(vn,'.');

K = 2*L+1;
plot(yn(4:end-3), '.'); % right now I am thinking that in order to get accurate symbol plots
% I go over the two ends the length of the filter. This is related to the
% ISI.
