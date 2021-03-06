% Template BER simulation
clearvars;
SIM = 'DFE';
% TODO: Now I would produce a plot comparing the trade off between
% performance and filter length. Loop over several curves.

SNR = 6:12; %dB
models = {'Ungerboeck', 'Forney', 'PSMFOM'};
ber = zeros(length(models), length(SNR));

for model_idx = 1:length(models)
    %% Load in an observation model -- Options ['Ungerboeck', 'Forney', 'PSMFOM']
    model = char(models(model_idx));

    load(['../Models/', model, '.mat']);
    L1 = eval([model, '.L1']);
    L2 = eval([model, '.L2']);
    hn = eval([model, '.filter']);
    if strcmp(model, 'Ungerboeck')
        noise_filter = eval([model, '.noise_filter']);
    end
    eval(['clear ' model]);
    fprintf(1, 'Simulating %s model\n', model);
    %%

    % Simulation data and parameters
    K1 = 5;
    K2 = L2; % We get no benefit out of making K1 larger than L2.
    filter_length = K1+K2+1;
    
    Nsymbols = 10e3;
    A = 1;

    LUT = A*[-1-1j;...
             -1+1j;...
             1-1j;...
             1+1j];

    %SNR = 5:12; %dB
    %ber = zeros(size(SNR));
    iters = 100;

    % Constants and Matricies that dont change for each iter.
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

    PHI = phi_first-phi_second;
    
    %%
    for idx = 1:length(SNR)

        EbNo = 10.^(SNR(idx)./10); % snr natural number
        var_n = A^2/(2*EbNo) * sum(abs(hn).^2) ;
        std_n = sqrt(var_n);
        
        err_cnt = 0;
        total_bits = 0;

        for idx2 = 1:iters
            % Generate Data
            symbol_idx = randi([1,4],Nsymbols,1);
            I = LUT(symbol_idx);

            s = conv(hn, I);
            s = s(-L1+1:end-L2);

            n = std_n*(randn(length(s),1) + 1j*randn(length(s),1));
            % create correlated noise for Ungerboeck
            if strcmp(model, 'Ungerboeck')
                n = filter(.83*noise_filter, 1, n);
            end

            vn = s + n;

            % Decision Feedback equalizer
            
            % Generate Noise covariance matrix
            N = convm(n, K1+1);
            Rn = zeros(K1+1, K1+1);
            for k = 1:K1+1
                for i = 1:K1+1
                    Rn(i, k) = N(:,i)'*N(:,k)/(Nsymbols);
                end
            end

            %  Solve wiener-hopf for feedforward
            Gamma = PHI + Rn;

            zeta = conj([zeros(K1-L2,1); flip(hn)]);
            zeta = zeta(1:K1+1);

            c_ff = Gamma\zeta;

            % solve for feedback
            hn_tmp = [zeros(K1+1-L2,1); flip(hn)];
            hn_tmp = hn_tmp(1:K1+1);

            H = convm(hn_tmp, L2).';
            H = H(1:L2, 1:K1+1);

            c_fb = -H*c_ff;
            
            % apply filter
            
            vn_tmp = [vn; zeros(K1,1)]; % append vn with zeros for feedforward terms.
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
            
            yn = I_tilde(L2+1:end);

            % detection
            err_cnt = err_cnt + sum(sign(real(yn)) ~= real(I));
            total_bits = total_bits + length(I);

            err_cnt = err_cnt + sum(sign(imag(yn)) ~= imag(I));
            total_bits = total_bits+ length(I);
        end
        ber(model_idx, idx) = err_cnt/total_bits;
        fprintf(1, 'SNR = %f dB, BER = %e\n', SNR(idx), ber(model_idx, idx));
    end
    % end models loop
end
%% Save results
SAVE = 0;
if SAVE
    save([SIM, '_BER.mat'], 'ber');
end

%% Plot
tX = 0:11;
tx = 10.^(0.1*tX);
ty = myQ(sqrt(2*tx));
figure(1); clf;
semilogy(SNR,ber,'o-',tX,ty,'--'); grid on;

title([SIM, ' Equalizer']);
xlabel('E_b/N_0 (dB)');
ylabel('BER');
legend('Ungerboeck (sim.)','Forney (sim.)', 'PSMFOM (sim.)','AWGN (theory)');

%% Print
printPlots = 0;
PATH = '/Users/mitchellburnett/Dropbox/BYU/776/Homework/final/img/';

if printPlots
    fig = figure(1);
    set(gca, 'FontSize', 12);
    width = 560; height = 160;
    xPos = 470; yPos = 980;
    set(fig, 'Position', [xPos, yPos, 600, 400]);
    fig.PaperPositionMode = 'auto';
    print(fig, [PATH, SIM, '_ber', '.eps'],'-depsc','-r0');
end
