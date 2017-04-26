function [ h, err ] = spike(g, n0, n)

% spike.m - Function for finding the FIR least squares inverse filter of 
%           g(n) for approximating a unit sample at time n = n0, i.e., 
%           delta(n-n0).
%           h - The coeff of the least squares approximation.
%           err - Approximation error.

    g = g(:);
    m = length(g);
    if (m+n-1) <= n0
        error('Delay too large')
    end
    
    G = convm(g,n);
    d = zeros(m+n-1,1);
    d(n0+1) = 1;
    h = G\d;
    err = 1 - G(n0+1,:)*h;
    
end

