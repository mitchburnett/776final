bits_per_symbol = 2;
LUT = [-1-1j, -1+1j, 1-1j, 1+1j];
M = 4;
L = 3;


% Gen States
Nstates = M^L;
STATES = zeros(Nstates,L);

dec_states = 0:Nstates-1;

binStates = dec2bin(dec_states);

for l = 1:Nstates
    bits = binStates(l,:);
    for m = 1:L
        idx = bin2dec(bits((m-1)*bits_per_symbol+1:m*bits_per_symbol));
        STATES(l,m) = LUT(idx+1);
    end
end

%% Init states -- the wind up on MLSE.
for k = 1:L-1
    n_init_states = M^k;
    %init_states = zeros(n_init_states,k); % maybe want L?
    init_states = zeros(n_init_states,k);
    
    dec_states = 0:n_init_states-1;
    
    binStates = dec2bin(dec_states);
    for l = 1:n_init_states
        bits = binStates(l,:);
        for m = 1:k
            idx = bin2dec(bits((m-1)*bits_per_symbol+1:m*bits_per_symbol));
            init_states(l,m) = LUT(idx+1);
        end
    end
    init_states;
%     state =  [zeros(n_init_states, L-k+1), init_states]; % prepend zeros
%     to have sizes correct for BM calculation.
%     x = [5/4, 1/2, 1/8, 1/32];
%     yk = -1/2;
%     BM = BM_Ungerboeck(state, x, -1/2);
    
end


%%

% Determine state row in matrix
cur_state = [STATES(14,:); STATES(1,:)];


%loop approach
L = length(cur_state);
decision = zeros(L,1);
for k = 1:L
    tmp = LUT - cur_state(k);
    tmp = tmp.^2;
    [~, decision(k)] = min(tmp);
end

% %matrix approach
% tmp = LUT-cur_state.';
% tmp = tmp.^2;
% [~, decision] = min(tmp,[],2);

% make decision zero based
decision = decision-1;

% create bit string
decode = dec2bin(decision,2);

% create a 1 dim array
%[xdim, ydim] = size(decode);
decode = reshape(decode', [1, length(decode)*bits_per_symbol]);

% Get the matlab based index of the state in the state matrix row. 
idx = bin2dec(decode) + 1;

