function [ states ] = gen_states(M, L, LUT, bits_per_symbol)

% I have found that states can be represented as a LUT in the sense that
% their symbol to bit mappings can represent the index that the state can
% be in a matric of states.

% This function creates the look up table for all possible states

    Nstates = M^L;
    
    states = zeros(Nstates, L);
    
    dec_states = 0:Nstates-1;
    
    bin_states = dec2bin(dec_states);
    
    for l = 1:Nstates
        bits = bin_states(l,:);
        for m = 1:L
            idx = bin2dec(bits((m-1)*bits_per_symbol+1:m*bits_per_symbol));
            states(l,m) = LUT(idx+1);
        end
    end

end

