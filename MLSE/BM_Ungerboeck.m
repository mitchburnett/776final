function [ BM ] = BM_Ungerboeck(state, x, yk)

% Branch metric for the Ungerboeck observation model in the viterbi 
% implementation of the MLSE.

%   There is a serial and vector mode for computation. I imagine the vector
%   mode will run faster in the end... if it runs at all...
    
     L = length(x) - 1;
     x = x(:); % make sure x is a column.
     
%     Ik = state(end);
%     Ik_prev = state(end-1:-1:end-L);
%     BM = Ik*(2*yk - Ik*x(1) - 2*Ik_prev*x(2:end));

    % how about vector operations....
    Ik = state(:,end);
    Ik_prev = state(:, end-1:-1:end-L);
    
    BM = real(conj(Ik).*(2*yk - Ik*x(1) - 2*Ik_prev*x(2:end)));
end

