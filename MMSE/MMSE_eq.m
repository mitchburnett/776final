% MMSE Equalizer
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

% Simulation data and parameters.
Nsymbols = 10e3;
A = 1;

LUT = A*[-1-1j;...
         -1+1j;...
         1-1j;...
         1+1j];
  
idx = randi([1,4],Nsymbols,1);
I = LUT(idx); % data

SNR = 20; %dB
EbNo = 10.^(SNR./10); % snr natural number
var_n = A^2/(2*EbNo) * sum(abs(hn).^2) ;
std_n = sqrt(var_n);


s = conv(hn, I);
%s = s(1:end-L2+1); % I still have a question as to why we dont cut off the front. Since a real data signal would not have come yet.
%s = s(-L1:end-L2); % And why here I didnt want to cut off anywhere but in zero forcing we wanted to??
n = std_n*(randn(length(s),1) + 1j*randn(length(s),1));

% create correlated noise for Ungerboeck
if strcmp(model, 'Ungerboeck')
    n = filter(.83*noise_filter, 1, n);
end

vn = s + n;

%% MMSE equalizer.
K1 = L2; % Get no benefit out of making K1 larger than L2. % Interesting... this is related to where the delta shows up in hn_eq = conv(w_mmse, hn); This may help me get closer to the answer to my question about how to pick out the data from the convolutions...
K2 = 30;
filter_length = K1+K2+1;

% Gamma matrix
H = convm(hn,filter_length);
Gamma = zeros(filter_length, filter_length);
for l = 1:filter_length
    for i = 1:filter_length
        Gamma(l,i) = H(:,i).'*conj(H(:,l)); % FN(:,i).'*conj(FN(:,l)) this seems to be the correct conjugate form.
    end
end

% Noise matrix
N = convm(n, filter_length);
Rn = zeros(filter_length, filter_length);
for l = 1:filter_length
    for i = 1:filter_length
        Rn(l,i) = N(:,i)'*N(:,l)/(Nsymbols);
    end
end

Gamma = Gamma + Rn;
zeta = conj([zeros(K1-L2,1); flip(hn); zeros(K2 + L1,1)]);

w_mmse = Gamma\zeta;

% apply filter
yn = conv(w_mmse, vn);
yn = yn(K1-L1+1:end-filter_length+1);

% equalized observation model.
hn_eq = conv(w_mmse, hn);
%%
nfft = 1024;
H_f = fftshift(fft(hn, nfft));
H_f_MMSE = fftshift(fft(w_mmse, nfft));
fspan = linspace(-.5, .5, nfft);

figure(1); clf;
plot(vn,'.');  grid on; hold on;
plot(yn, '.');

figure(2); clf;
stem(0:length(hn_eq)-1, abs(hn_eq)); grid on;

figure(3); clf;
plot(fspan, 20*log10(abs(H_f))); grid on; hold on;
plot(fspan, 20*log10(abs(H_f_MMSE)));
legend('F', 'MMSE');

