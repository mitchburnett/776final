% DF Equalizer Code
clearvars;

%% Load in an observation model -- Options ['Ungerboeck', 'Forney', 'PSMFOM']
model = 'PSMFOM';

load(['../Models/', model, '.mat']);
L1 = eval([model, '.L1']);
L2 = eval([model, '.L2']);
hn = eval([model, '.filter']);
eval(['clear ' model]);
% set up symbols
Nsymbols = 10e3;
%Nsymbols = 100;
A = 1;

LUT = A*[-1-1j;...
         -1+1j;...
         1-1j;...
         1+1j];
  
idx = randi([1,4],Nsymbols,1);
I = LUT(idx);

SNR = 55; %dB
EbNo = 10.^(SNR./10); %SNR
sigma_n = sqrt(A^2/(2*EbNo));

%n = sigma_n*(randn(Nsymbols,1) + 1j*randn(Nsymbols,1));
%r = I + n;

IF = conv(hn, I);
%IF = IF(1:end-(length(hn)-1));
%IF = IF(1:end-L2);
%IF = IF(-L1+1:end);

% Forney model output
n = sigma_n*(randn(length(IF),1) + 1j*randn(length(IF),1));
vn = IF + n;

%% Setup to solve for filter coeff

% Total filter length is K1+L2+1. DFE filter has support -K1 <= n <= L2
K1 = 5;
K2 = L2; % We get no benefit out of making K1 larger than L2.
filter_length = K1+K2+1;

% Generate phi matrix.
H = convm(hn,K1+1);
Hflip = flip(H, 2); %flip the matrix around vertical center.

phi_first = zeros(K1+1, K1+1);
for k = -K1:0
    for i = -K1:0
        row_start = k+K1+1;
        row_end = row_start+(L2-L1);
        % pulls the column for the first conj term in phi.
        tmp = Hflip(row_start:row_end,-i+1);
        phi_first(K1+i+1,K1+k+1) = hn.'*conj(tmp);
    end
end

% same goes for phi_second. I didnt even try using convm. However, it may
% be worth seeing if there is a more tractable way.
phi_second = zeros(K1+1, K1+1);
hn_tmp = [hn; zeros(K1, 1)];
for k = -K1:0
    for i = -K1:0
        % notice the +2. There is an extra +1 since the second sum in phi
        % starts from l=1 and extends to L2.
        h1 = hn_tmp(-k+2-L1:-k+L2+1-L1);
        h2 = hn_tmp(-i+2-L1:-i+L2+1-L1);
        phi_second(K1+i+1,K1+k+1) = h1.'*conj(h2);
    end
end

phi = phi_first-phi_second;

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

zeta = conj([zeros(K1-L2,1); flip(hn)]);
zeta = zeta(1:K1+1);

c_ff = phi\zeta;

% Feedback coeff
hn_tmp = [zeros(K1+1-L2,1); flip(hn)];
hn_tmp = hn_tmp(1:K1+1);

H = convm(hn_tmp, L2).';
H = H(1:L2, 1:K1+1);

c_fb = -H*c_ff;

%% Plot channel and equalizer response

ck = [c_ff; c_fb];

nfft = 1024;
H_f = fftshift(fft(hn, nfft));
H_F_DF = fftshift(fft(ck, nfft));
fspan = linspace(-.5, .5, nfft);

figure(3); clf;
plot(fspan, 20*log10(abs(H_f))); grid on; hold on;
plot(fspan, 20*log10(abs(H_F_DF)));
legend('hn', 'DFE');
title('Channel Response');

hn_eq = conv(ck, hn);
figure(5); clf;
stem(0:length(hn_eq)-1, abs(hn_eq)); grid on;

%% Apply Filter

% append vn with zeros for feedforward terms.
vn_tmp = conv(vn,c_ff);
I_tilde = zeros(Nsymbols + L2, 1);
I_out = zeros(Nsymbols,1);
I_hat = zeros(Nsymbols,1);

for n = 1:Nsymbols
    % Calculate I_hat
    %vn_vec = vn_tmp(n:K1+n);
    vn_vec = vn_tmp(n);
    I_tilde_vec = I_tilde(n:n+L2-1);
    %I_hat(n) = c_ff.'*vn_vec + c_fb.'*I_tilde_vec(L2:-1:1);
    I_hat(n) = vn_vec + c_fb.'*I_tilde_vec(L2:-1:1);
    
    % make decision for I_tilde
    tmp_d = LUT - I_hat(n);
    tmp_d = tmp_d.^2;
    [~, decision] = min(tmp_d);
    
    I_out(n) = decision;
    I_tilde(n+L2) = LUT(decision);
end

BER = sum(I_out==idx)/Nsymbols;

%% Plot constellation
figure(1); clf;
plot(vn,'.');  grid on; hold on;
plot(I_hat, '.');
legend('v_n', 'I_{hat}');
axis('square');






