% MLSE BER simulation
clearvars;
SIM = 'MLSE';

%models = {'Ungerboeck', 'Forney', 'PSMFOM'};
models = {'Ungerboeck', 'Forney'}; % add PSMFOM back into the legend.

% Global simulation variables and parameter constants
SNR = 6:9; %dB
ber = zeros(length(models), length(SNR));

Nsymbols = 10e1;
A = 1;

LUT = A*[-1-1j;...
         -1+1j;...
         1-1j;...
         1+1j];
    
bits_per_symbol = 2;
M = 4;
iters = 25;

for model_idx = 1:length(models)
    %% Load in an observation model -- Options ['Ungerboeck', 'Forney', 'PSMFOM']
    model = char(models(model_idx));

    load(['../Models/', model, '.mat']);
    L1 = eval([model, '.L1']);
    L2 = eval([model, '.L2']);
    hn = eval([model, '.filter']);
    if strcmp(model, 'Ungerboeck')
        metric = 1;
        x = hn;
        hn = hn(-L1+1:end);
    
        L = L2;
        noise_filter = eval([model, '.noise_filter']);
    else
        L = -L1+L2;
    end
    eval(['clear ' model]);
    fprintf(1, 'Simulating %s model\n', model);
    
    %%
    % Generate a matrix of all possible states
    Nstates = M^L;
    STATES = gen_states(M, L, LUT, bits_per_symbol); % (NStates x L) matrix of possible states.

    for idx = 1:length(SNR)

        EbNo = 10.^(SNR(idx)./10); % snr natural number
        var_n = A^2/(2*EbNo) * sum(abs(hn).^2) ;
        std_n = sqrt(var_n);
        
        err_cnt = 0;
        total_bits = 0;

        for idx2 = 1:iters
            fprintf(1, 'iter = %d/%d\n', idx2, iters);
            % Generate info and winner matrix for current iteration
            INFO = zeros(Nstates, Nsymbols);%, 2); % 3-d matrix. First z is for them PM. Second z is for the idx of the winner?
            WINNER = zeros(Nstates, Nsymbols);
            
            % Generate Data
            symbol_idx = randi([1,M],Nsymbols,1);
            I = LUT(symbol_idx);
            
            if strcmp(model, 'Ungerboeck') || strcmp(model, 'Dummy')
                s = conv(I, x);
            else
                s = conv(I, hn);
            end
            s = s(-L1+1:end-L2);

            n = std_n*(randn(length(s),1) + 1j*randn(length(s),1));
            % create correlated noise for Ungerboeck
            if strcmp(model, 'Ungerboeck')
                n = filter(.83*noise_filter, 1, n);
            end

            vn = s + n;

            % MLSE Equalizer
            
            % Init the algorithm -- i.e, the wind up.
            init_PM = zeros(1,1);
            for k = 1:L
                init_states = gen_states(M,k, LUT,bits_per_symbol);
                n_init_states = M^k;

                % Prepend with zeros to match branch metric calculation.
                state = [zeros(n_init_states, L-k+1), init_states];

                if strcmp(model, 'Ungerboeck')
                    BM = BM_Ungerboeck(state, hn, vn(k));
                elseif strcmp(model, 'Forney')
                    BM = -BM_Forney(state, hn, vn(k));
                elseif strcmp(model, 'PSMFOM')
                    BM = -BM_PSMFOM(state, hn, vn(k), L1);
                end

                init_PM = repmat(init_PM, [1, M]);
                init_PM = reshape(init_PM', [n_init_states,1]);
                init_PM = init_PM + BM;
            end

            INFO(:, k) = init_PM;
            
            % and we are off...
            
            full_state = zeros(M, L+1);
            prev_state_idx = zeros(M,1);
            for k = L+1:Nsymbols
                for l = 1:Nstates
                    cur_state = STATES(l,:);

                    full_state(:, 2:end) = repmat(cur_state, [M,1]);
                    full_state(:, 1) = LUT;

                    for m = 1:M
                        prev_state_idx(m) = get_state_idx(full_state(m, 1:L), LUT, bits_per_symbol);
                    end

                    tmpPM = INFO(prev_state_idx, k-1);

                    if strcmp(model, 'Ungerboeck')
                        BM = BM_Ungerboeck(full_state, hn, vn(k));
                    elseif strcmp(model, 'Forney')
                        BM = -BM_Forney(full_state, hn, vn(k));
                    elseif strcmp(model, 'PSMFOM')
                        BM = -BM_PSMFOM(full_state, hn, vn(k), L1);
                    end

                    tmpPM = tmpPM + BM;

                    [s, sort_idx] = sort(tmpPM, 'descend');

                    winner_state_idx = get_state_idx(full_state(sort_idx(1), 1:L), LUT, bits_per_symbol);

                    INFO(l, k) = tmpPM(sort_idx(1));
                    WINNER(l, k) = winner_state_idx;
                end
            end
            
            % let's Unwind...
            [X, UIDX] = sort(INFO, 'descend');

            max_idx = UIDX(1,:);
            max_idx = max_idx(L:end);

            I_hat = zeros(length(max_idx),1);
            for n = length(max_idx):-1:1
                state = STATES(max_idx(n),:);

                decision = state(end);

                I_hat(n) = decision;

            end
            
            % unwind based on 'who got me to here'
            I_hat2 = zeros(length(max_idx),1);
            state = STATES(max_idx(end),:);
            decision = state(end);
            I_hat2(end) = decision;
            who_got_me_here = max_idx(end);
            for n = length(max_idx):-1:2
                who_got_me_here = WINNER(who_got_me_here, n+ (L-1)); % L-1 because we chopped some k's off in the beginning and I shorten max_idx by this.
                state = STATES(who_got_me_here, :);
                decision = state(end);
                I_hat2(n-1) = decision;
            end

            I = I(L:end);
            
            % detection
            err_cnt = err_cnt + sum(sign(real(I_hat)) ~= real(I));
            total_bits = total_bits + length(I);

            err_cnt = err_cnt + sum(sign(imag(I_hat)) ~= imag(I));
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
%legend('Ungerboeck (sim.)','Forney (sim.)', 'PSMFOM (sim.)','AWGN (theory)');
legend('Ungerboeck (sim.)','Forney (sim.)','AWGN (theory)');

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
