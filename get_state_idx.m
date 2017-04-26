function [ idx ] = get_state_idx(cur_state, LUT, bits_per_symbol)

% This function extend the state to bit LUT method, to where given a state,
% this function determins the index in the state matrix that it corresponds
% to.

%   loop approach
    L = length(cur_state);
    decision = zeros(L,1);
    for k = 1:L
        tmp = LUT - cur_state(k);
        tmp = tmp.^2;
        [~, decision(k)] = min(tmp);
    end

%     %matrix approach matlab 2016
%     tmp = LUT-cur_state.';
%     tmp = tmp.^2;
%     [~, decision] = min(tmp,[],2);

    % make decision zero based
    decision = decision-1;

    % create bit string
    decode = dec2bin(decision,bits_per_symbol);

    % create a 1 dim array
    %[xdim, ydim] = size(decode);
    decode = reshape(decode', [1, length(decode)*bits_per_symbol]);

    % Get the matlab based index of the state in the state matrix row. 
    idx = bin2dec(decode) + 1;


end

