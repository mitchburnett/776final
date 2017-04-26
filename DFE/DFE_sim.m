% Decision Feedback Equalizer BER simulation
clearvars;

% TODO: Now I would produce a plot comparing the trade off between
% performance and filter length. Loop over several curves.

%% Load in an observation model -- Options ['Ungerboeck', 'Forney', 'PSMFOM']
model = 'PSMFOM';

load(['../Models/', model, '.mat']);
L1 = eval([model, '.L1']);
L2 = eval([model, '.L2']);
hn = eval([model, '.filter']);
eval(['clear ' model]);
%%

% Simulation data and parameters
K1 = L2; % We get no benefit out of making K1 larger than L2.
K2 = 5;
filter_length = K1+K2+1;

Nsymbols = 10e4;
A = 1;

LUT = A*[-1-1j;...
         -1+1j;...
         1-1j;...
         1+1j];
     
SNR = 5:12; %dB
ber = zeros(size(SNR));
iters = 100;

% Gamma matrix - doesnt need to be updated, unlike noise matrix.
H = convm(hn,filter_length);
Gamma = zeros(filter_length, filter_length);
for l = 1:filter_length
    for i = 1:filter_length
        Gamma(l,i) = H(:,i).'*conj(H(:,l));
    end
end
%%
for idx = 1:length(SNR)
    
    EbNo = 10.^(SNR(idx)./10); % snr natural number
    sigma_n = sqrt(A^2/(2*EbNo));
    
    err_cnt = 0;
    total_bits = 0;
    
    for idx2 = 1:iters
        symbol_idx = randi([1,4],Nsymbols,1);
        I = LUT(symbol_idx);

        s = conv(hn, I);
        %s = s(1:end-L2+1);
        
        n = sigma_n*(randn(length(s),1) + 1j*randn(length(s),1));
        
        vn = s + n;

        % MMSE Equalizer

        % Generate noise covariance matrix
        N = convm(n, filter_length);
        Rn = zeros(filter_length, filter_length);
        for l = 1:filter_length
            for i = 1:filter_length
                Rn(l,i) = N(:,i)'*N(:,l)/(Nsymbols);
            end
        end

        % Solve wiener-hopf
        Gamma = Gamma + Rn;
        zeta = conj([zeros(K1-L2,1); flip(hn); zeros(K2 + L1,1)]);

        w_mmse = Gamma\zeta;

        % apply filter
        yn = conv(w_mmse, vn);
        yn = yn(K1-L1+1:end-filter_length+1);
        
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
