% DF Equalizer Code - testing Forney observation model.

clearvars;
% load in an observation model.
load('../Models/Forney.mat');
fn = Forney.filter.'; clear Forney;
L = length(fn)-1;

% set up symbols
Nsymbols = 10e3;
%Nsymbols = 30;
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

%n = sigma_n*(randn(Nsymbols,1) + 1j*randn(Nsymbols,1));
%r = I + n;

IF = conv(fn, I);
IF = IF(1:end-(length(fn)-1));

% Forney model output
n = sigma_n*(randn(length(IF),1) + 1j*randn(length(IF),1));
vn = IF + n;

%% Setup to solve for filter coeff

% Total filter length is K1+L2+1
K1 = 10;

% General observation model. Forney model L1 = 0, L2 = L.
L1 = 0;
L2 = L;

% % Generate phi matrix.
% H = convm(fn,K1+1);
% phi1 = zeros(K1+1, K1+1);
% for k = 1:K1+1
%     for i = 1:K1+1
%         phi1(i,k) = H(:,i)'*H(:,k);
%     end
% end
% 
% H2 = convm(fn(2:end), K1+1);
% phi2 = zeros(K1+1, K1+1);
% for k = 1:K1+1
%     for i = 1:K1+1
%         phi2(i,k) = H2(L2+1:end,i)'*H2(L2+1:end,k);
%     end
% end
% 
% phi = phi1 - phi2;
% 

% fn_tmp = [zeros(1, K1), fn, zeros(1,K1-L)]; % pad fn on both sides since the sum goes to K1 this will cause it to go beyond its support.
% phi_forney = zeros(K1+1, K1+1);
% for k = -K1:0
%     for i = -K1:0
%         % this formulation the sum reaches negative indices and so the
%         % added K1+1 in tmp1 and tmp2 are to center back at 0. The K1+1 in
%         % the phi_forney assignment is to have the top left of phi(-K1,-K1)
%         % map to (1,1).
%         tmp1 = fn_tmp(K1+1:-i+K1+1);
%         tmp2 = fn_tmp(i-k+K1+1:i-k+K1+1-i);
%         phi_forney(K1+i+1,K1+k+1) = conj(tmp1)*tmp2.';
%     end
% end
% phi = phi_forney;

H = convm(fn,K1+1);
Hflip = flip(H, 2); %flip the matrix around vertical center.

phi_first = zeros(K1+1, K1+1);
for k = -K1:0
    for i = -K1:0
        row_start = k+K1+1;
        row_end = row_start+L2;
        % pulls the column for the first conj term in phi.
        tmp = Hflip(row_start:row_end,-i+1);
        phi_first(K1+i+1,K1+k+1) = fn*conj(tmp);
    end
end

% same goes for phi_second. I didnt even try using convm. However, it may
% be worth seeing if there is a more tractable way.
phi_second = zeros(K1+1, K1+1);
fn_tmp = [fn, zeros(1, K1)];
for k = -K1:0
    for i = -K1:0
        % notice the +2. There is an extra +1 since the second sum in phi
        % starts from l=1 and extends to L2.
        h1 = fn_tmp(-k+2:-k+L2+1);
        h2 = fn_tmp(-i+2:-i+L2+1);
        phi_second(K1+i+1,K1+k+1) = h1*h2';
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

zeta = [zeros(1,K1-L), flip(fn)]';

c_ff = phi\zeta;

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

figure(3); clf;
plot(fspan, 20*log10(abs(F))); grid on; hold on;
plot(fspan, 20*log10(abs(F_DF)));
legend('Fn', 'DFE');
title('Channel Response');

%% Apply Filter

% append vn with zeros for feedforward terms.
vn_tmp = [vn; zeros(K1,1)];
c_ff = flip(c_ff); % flip c_ff for conveience such that rvec increases forward.

I_tilde = zeros(Nsymbols + L2, 1);
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
    I_tilde(n+L2) = LUT(decision);
end

BER = sum(I_out==idx)/Nsymbols;

%% Plot constellation
figure(1); clf;
plot(vn,'.');  grid on; hold on;
plot(I_hat, '.');
legend('v_n', 'I_{hat}');
axis('square');

hn_eq = conv(ck, fn);
figure(5); clf;
stem(0:length(hn_eq)-1, abs(hn_eq)); grid on;






