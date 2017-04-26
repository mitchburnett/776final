clearvars;

%%Load in an observation model -- Options ['Ungerboeck', 'Forney', 'PSMFOM', 'Dummy']
model = 'Ungerboeck';

load(['../Models/', model, '.mat']);
L1 = eval([model, '.L1']);
L2 = eval([model, '.L2']);
hn = eval([model, '.filter']);
% set the branch metric to use
if strcmp(model, 'Ungerboeck') || strcmp(model, 'Dummy')
    x = hn;
    hn = hn(-L1+1:end);
    
    L = L2;
else
    L = -L1+L2;
end
eval(['clear ' model]);

%%
% TEST PARAMS
NOISE = 1;

% Parameters
A = 1;
LUT = A*[-1-1j; -1+1j; 1-1j; 1+1j];
%LUT = A*[-1; 1];
bits_per_symbol = 2;
M = 4;
Nstates = M^L;

% Generate data
Nsymbols = 1000;

SNR = 9; %dB
EbNo = 10.^(SNR./10); % snr natural number

idx = randi([1,M],Nsymbols,1);
I = LUT(idx);

if strcmp(model, 'Ungerboeck') || strcmp(model, 'Dummy')
    yn = conv(I, x);
else
    yn = conv(I, hn);
end

if NOISE
    sigma_n = sqrt(A^2/(2*EbNo));
else
    sigma_n = 0;
end

yn = yn(-L1+1:end-L2);
n = sigma_n*sigma_n*(randn(length(yn),1) + 1j*randn(length(yn),1));

yn = yn + n;

%%

% Generate a matrix of all possible states
STATES = gen_states(M, L, LUT, bits_per_symbol); % (NStates x L) matrix of possible states.
INFO = zeros(Nstates, Nsymbols);%, 2); % 3-d matrix. First z is for them PM. Second z is for the idx of the winner?
WINNER = zeros(Nstates, Nsymbols);
% Init the algorithm -- i.e, the wind up.
init_PM = zeros(1,1);
for k = 1:L
    init_states = gen_states(M,k, LUT,bits_per_symbol);
    n_init_states = M^k;
    
    % Prepend with zeros to match branch metric calculation.
    state = [zeros(n_init_states, L-k+1), init_states];
    
    if strcmp(model, 'Ungerboeck')
        BM = BM_Ungerboeck(state, hn, yn(k));
    elseif strcmp(model, 'Forney')
        BM = -BM_Forney(state, hn, yn(k));
    elseif strcmp(model, 'PSMFOM')
        BM = -BM_PSMFOM(state, hn, yn(k), L1);
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
            BM = BM_Ungerboeck(full_state, hn, yn(k));
        elseif strcmp(model, 'Forney')
            BM = -BM_Forney(full_state, hn, yn(k));
        elseif strcmp(model, 'PSMFOM')
            BM = -BM_PSMFOM(full_state, hn, yn(k), L1);
        end

        tmpPM = tmpPM + BM;

        [s, idx] = sort(tmpPM, 'descend');

        winner_state_idx = get_state_idx(full_state(idx(1), 1:L), LUT, bits_per_symbol);
        
        INFO(l, k) = tmpPM(idx(1));
        WINNER(l, k) = winner_state_idx;
    end
end

%% lets Unwind...
[X, IDX] = sort(INFO, 'descend');

max_idx = IDX(1,:);
max_idx = max_idx(L:end);

I_hat = zeros(length(max_idx),1);
for n = length(max_idx):-1:1
    state = STATES(max_idx(n),:);
    
    decision = state(end);
    
    I_hat(n) = decision;
    
end

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

err_real = sum(sign(real(I_hat)) ~= real(I))
err_imag = sum(sign(imag(I_hat)) ~= imag(I))














