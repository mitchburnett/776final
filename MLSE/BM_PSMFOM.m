function [ BM ] = BM_PSMFOM(state, h, vk, L1)

% Branch metri for the Forney observation model in the viterbi
% implementation of the MLSE.
    
    L = length(h) - 1;
    h = h(:);
    Ik = state(:, end:-1:end-L);
    
    Ik = circshift(Ik, -L1, 2); % shift to the right to line up.
    BM = abs(vk - Ik*h).^2;


end

