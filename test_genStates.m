bits_per_symbol = 2;
LUT = [-1-1j, -1+1j, 1-1j, 1+1j];
M = 4;
L = 3;

% Gen States
STATES = gen_states(M, L, LUT, bits_per_symbol);

%% Init states -- the windup on MLSE.
for k = 1:L-1
    init_state = gen_states(M,k, LUT,bits_per_symbol);
end

%% Determine state row in matrix
cur_state = STATES(6,:);
idx = get_state_idx(cur_state, LUT, bits_per_symbol);