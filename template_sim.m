% Template BER simulation
clearvars;
SIM = 'DFE';
% TODO: Now I would produce a plot comparing the trade off between
% performance and filter length. Loop over several curves.

SNR = 6:18; %dB
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

    Nsymbols = 10e3;
    A = 1;

    LUT = A*[-1-1j;...
             -1+1j;...
             1-1j;...
             1+1j];

    %SNR = 5:12; %dB
    %ber = zeros(size(SNR));
    iters = 100;

    % Constants and Matricies that dont change for each iter, but are
    % specific to the model.
    
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
            %s = s(1:end-L2+1);

            n = std_n*(randn(length(s),1) + 1j*randn(length(s),1));
            % create correlated noise for Ungerboeck
            if strcmp(model, 'Ungerboeck')
                n = filter(.83*noise_filter, 1, n);
            end

            vn = s + n;

            % Equalizer


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
SAVE = 1;
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
printPlots = 1;
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
