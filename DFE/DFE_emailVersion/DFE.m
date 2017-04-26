% DF Equalizer Code - testing Forney observation model.

clearvars;
% load in an observation model.
load('ForneyModel.mat');
fn = f_n; clear f_n;
L = length(fn)-1;

% set up symbols
Nsymbols = 10e2;
A = 1;

LUT = A*[-1-1j;...
         -1+1j;...
         1-1j;...
         1+1j];
  
idx = randi([1,4],Nsymbols,1);
I = LUT(idx);

SNR = 25; %dB
EbNo = 10.^(SNR./10); %SNR
sigma_n = sqrt(A^2/(2*EbNo));

n = sigma_n*(randn(Nsymbols,1) + 1j*randn(Nsymbols,1));
r = I + n;

IF = conv(fn, I);
IF = IF(1:end-(length(fn)-1));

% Forney model output
vn = IF + n;

%% Setup to solve for filter coeff

% Total filter length is K1+L2+1
K1 = 3;

% General observation model. Forney model L1 = 0, L2 = L.
L1 = 0;
L2 = L;

% Generate phi matrix.
H = convm(fn,K1+1);
phi1 = zeros(K1+1, K1+1);
for k = 1:K1+1
    for i = 1:K1+1
        phi1(i,k) = H(:,i)'*H(:,k);
    end
end

H2 = convm(fn(2:end), K1+1);
phi2 = zeros(K1+1, K1+1);
for k = 1:K1+1
    for i = 1:K1+1
        phi2(i,k) = H2(L2+1:end,i)'*H2(L2+1:end,k);
    end
end

phi = phi1 - phi2;

% Generate Noise covariance matrix
N = convm(n, K1+1);
Rn = zeros(K1+1, K1+1);
for k = 1:K1+1
    for i = 1:K1+1
        Rn(i, k) = N(:,i)'*N(:,k)/(Nsymbols);
    end
end

%% Calculate coefficients

% Feed-forward coeff
phi = phi + Rn;

h_vec = [zeros(1,K1-L), flip(fn)]';

c_ff = phi\h_vec;

% Feedback coeff
fn_tmp = [zeros(1, K1+1-L), flip(fn)];
fn_tmp = fn_tmp(1:K1+1);

H = convm(fn_tmp, L2).';
H = H(1:L2, 1:K1+1);

c_fb = -H*c_ff;

%% Plot channel and equalizer response

ck = [c_ff; c_fb];

nfft = 1024;
F = fftshift(fft(fn, nfft));
F_DF = fftshift(fft(ck, nfft));
fspan = linspace(-.5, .5, nfft);

figure(1); clf;
plot(fspan, 20*log10(abs(F))); grid on; hold on;
plot(fspan, 20*log10(abs(F_DF)));
legend('Fn', 'DFE');
title('Channel Vs. Equalizer Response');
xlabel('Normalized Frequency (rad/sample)');
ylabel('Magnitude (dB)');

%% Apply Filter

% append vn with zeros for feedforward terms.
vn_tmp = [vn; zeros(K1,1)];
c_ff = flip(c_ff); % flip c_ff for conveience such that rvec increases forward.

I_tilde = zeros(Nsymbols + K1, 1);
I_out = zeros(Nsymbols,1);
I_hat = zeros(Nsymbols,1);

for n = 1:Nsymbols
    % Calculate I_hat
    vn_vec = vn_tmp(n:K1+n);
    I_tilde_vec = I_tilde(n:n+L2-1);
    I_hat(n) = c_ff.'*vn_vec + c_fb.'*I_tilde_vec(L2:-1:1);
    
    % make decision for I_tilde
    tmp_d = LUT - I_hat(n);
    tmp_d = tmp_d.^2;
    [~, decision] = min(tmp_d);
    
    I_out(n) = decision;
    I_tilde(n+K1) = LUT(decision);
end

BER = sum(I_out==idx)/Nsymbols;

%% Plot constellation
figure(2); clf;
plot(vn,'.');  grid on; hold on;
plot(I_hat, '.');
legend('v_n', 'I_{hat}');
axis('square');
xlabel('Real');
ylabel('Imaginary');

