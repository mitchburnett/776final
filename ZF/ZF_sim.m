% Zero Forcing Equalizer BER Simulation
clearvars;

% TODO: Now I would produce a plot comparing the trade off between
% performance and filter length. Loop over several curves.

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
P = 40;  % Equalizer length.
n0 = 15; % position of delta in equalized signal. Needs to be >= -L1.
Nsymbols = 10e4;
A = 1;

LUT = A*[-1-1j;...
         -1+1j;...
         1-1j;...
         1+1j];
     
SNR = 5:18; %dB
ber = zeros(size(SNR));
iters = 100;   

for idx = 1:length(SNR)
    
    EbNo = 10.^(SNR(idx)./10); % snr natural number
    var_n = A^2/(2*EbNo) * sum(abs(hn).^2) ;
    std_n = sqrt(var_n);
    
    err_cnt = 0;
    total_bits = 0;
    
    for idx2 = 1:iters
        symbol_idx = randi([1,4],Nsymbols,1);
        I = LUT(symbol_idx);

        s = conv(hn, I);
        s = s(1:end-L2+1);
        
        n = std_n*(randn(length(s),1) + 1j*randn(length(s),1));
        
        % create correlated noise for Ungerboeck
        if strcmp(model, 'Ungerboeck')
            n = filter(.83*noise_filter, 1, n);
        end
        
        vn = s + n;

        % ZF equalizer
        [w_zf, err] = spike(hn, n0, P);
        
        % apply filter
        yn = conv(w_zf, vn);
        yn = yn(n0+1:end-P+n0+L1);
        
        % detection
        err_cnt = err_cnt + sum(sign(real(yn)) ~= real(I));
        total_bits = total_bits + length(I);
        
        err_cnt = err_cnt + sum(sign(imag(yn)) ~= imag(I));
        total_bits = total_bits+ length(I);
    end
    ber(idx) = err_cnt/total_bits;
    fprintf(1, 'SNR = %f dB, BER = %e\n', SNR(idx), ber(idx));
end

tX = 0:11;
tx = 10.^(0.1*tX);
ty = myQ(sqrt(2*tx));
figure(776); clf;
semilogy(SNR,ber,'o-',tX,ty,'--'); grid on;
xlabel('E_b/N_0 (dB)');
ylabel('BER');
legend('QPSK (sim.)','AWGN (theory)');
