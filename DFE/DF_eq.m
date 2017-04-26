% This code was my first attempt at the DFE.

% DF Equalizer

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
I = LUT(idx); % data

SNR = 25; %dB
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

%% DF Equalizer

% TODO: I still need to vet my algorithm here for computing
% the matrix. I knokw one issue is I am putting phi(-K1, -K1)
% in the bottom right entry of the matrix instead of of the first So I
% believe I need to have (i+K1+1, k+K1+1) in the assignment instead.

% FOLLOW UP: Gamma2 seems to be more correct.
% The above todo talks about using (i+K1+1, k+K1+1) and I have
% done a new script that computes the phi used in the DFE and this
% assignment seems to work better.

% ALSO: Although Gamma2 may be close to right, the original handout had an
% index error where the algorithm below wouldnt calculate tmp2 and tmp3
% correctly. tmp1 may be correct.

% solve for feedforward coeff.
K1 = 3;
FN = convm(fn,K1+1);

Gamma = zeros(K1+1, K1+1);
Gamma2 = zeros(K1+1, K1+1);

fn_ex = [zeros(1, K1), fn, zeros(1, K1-L)];
fn_ex2 = [fn(2:end), zeros(1,K1)];
for k = -K1:0
    for i = -K1:0
        tmp1 = FN(:,-i+1).'*conj(FN(:,-k+1));
        tmp2 = fn_ex((i-k)+K1+1)*ones(1,L);
        tmp3 = fn_ex2((-i+1):(-i+L));
        Gamma(-i+1,-k+1) = tmp1 - tmp2*tmp3';
        Gamma2(i+K1+1, k+K1+1) = tmp1 - tmp2*tmp3';
    end
end

keyboard;
%%
N = convm(n, K1+1);
Rn = zeros(K1+1, K1+1);
for k = 1:K1+1
    for i = 1:K1+1
        Rn(i, k) = N(:,i)'*N(:,k)/(Nsymbols);
    end
end

Gamma2 = Gamma2 + Rn;

zeta = [zeros(1,K1-L), flip(fn)]';

c_ff = Gamma2\zeta;

%% fb coeff

fn_ex = [fn, zeros(1, K1+1-L)];
fn_ex = flip(fn_ex);
fn_ex = fn_ex(1:K1+1);

H = convm(fn_ex, L).';
H = H(1:L, 1:K1+1);

c_fb = -H*c_ff;

%% apply feed forward

% FOLLOW UP: This calculation is not correct. I need to use an adaptice
% filter sort of implementation since this is a non linear filter.

yn = conv(c_ff, vn);

figure(1); clf;
%plot(yn(16:end-14), '.');
plot(I, '.'); grid on; hold on;
plot(r, '.');
plot(vn,'.');

plot(yn , '.');

% create decisions
In_hat = zeros(length(yn),1);
for k = 1:length(yn)
    tmp_decision = LUT - yn(k);
    tmp_decision = tmp_decision.^2;
    [~, In_hat(k)] = min(tmp_decision);
end

% make In_tilde signal points again.
In_hat = LUT(In_hat);

% apply feedback
yn_hat = conv(c_fb, In_hat);

I_tilde = yn + yn_hat(3:end);








