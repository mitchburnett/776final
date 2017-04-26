clearvars;

%%Load in an observation model -- Options ['Ungerboeck', 'Forney', 'PSMFOM', 'Dummy']
model = 'Forney';

load(['Models/', model, '.mat']);
FL1 = eval([model, '.L1']);
FL2 = eval([model, '.L2']);
f = eval([model, '.filter']);

model = 'Ungerboeck';
load(['Models/', model, '.mat']);
XL1 = eval([model, '.L1']);
XL2 = eval([model, '.L2']);
x = eval([model, '.filter']);

% Generate data
Nsymbols = 10e5;
M = 4;
A = 1;
LUT = A*[-1-1j; -1+1j; 1-1j; 1+1j];
idx = randi([1,M],Nsymbols,1);
I = LUT(idx);

SNR = 5; %dB
EbNo = 10.^(SNR./10); % snr natural number

x_power = sum(abs(x).^2);

newx = conv(f, conj(flip(f)));


x0 = real(x(-XL1+1));
%var_n = A^2/(2*EbNo) * sum(abs(newx).^2);% * x_power; % Even though it is real I took the real to cast to a double.
var_n = A^2/(2*EbNo);

std_n = sqrt(var_n);

v = conv(I, newx);
v = v(-XL1+1:end-XL2);

n = std_n*(randn(Nsymbols,1) + 1j*randn(Nsymbols,1));

n = filter(f, 1, n);

r = v + n;




