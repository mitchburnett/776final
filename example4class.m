%%% example4class

LUT = [-1-1i; -1+1i; 1-1i; 1+1i];
%SNR = 3:10;
SNR = 10:10;
ber = zeros(size(SNR));
its = 100;

for idx = 1:length(SNR)
    snr = 10^(0.1*SNR(idx));
    nvar = 0.5/snr;
    nstd = sqrt(nvar);
    err_cnt = 0;
    total_bits = 0;
    for idx2 = 1:its
        ii = randi([0 3],100000,1);
        I = LUT(ii+1);
        y = I + nstd*(randn(size(I))+1i*randn(size(I)));
        err_cnt = err_cnt + sum(sign(real(y)) ~= real(I));
        total_bits = total_bits + length(I);
        err_cnt = err_cnt + sum(sign(imag(y)) ~= imag(I));
        total_bits = total_bits + length(I);
    end
    ber(idx) = err_cnt/total_bits;
    fprintf(1,'SNR = %f dB, BER = %e\n',SNR(idx),ber(idx));
end

tX = 0:11;
tx = 10.^(0.1*tX);
ty = myQ(sqrt(2*tx));
figure(776); clf;
semilogy(SNR,ber,'o-',tX,ty,'--'); grid on;
xlabel('E_b/N_0 (dB)');
ylabel('BER');
legend('QPSK (sim.)','AWGN (theory)');
