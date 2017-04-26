function [ BM ] = BM_Forney(state, f, vk)

% Branch metri for the Forney observation model in the viterbi
% implementation of the MLSE.
    
    L = length(f) - 1;
    f = f(:);
    Ik = state(:, end:-1:end-L);
    BM = abs(vk - Ik*f).^2;


end

