% experimenting with how to build the convolution matrix for DF.

clearvars;
L1 = 1;
L2 = 2;
K1 = 3;

    %  -1  0   1    2
h = [-1/2, 1, 1/2, 1/3];
h_ex = [h, zeros(1,K1+1-L2)];
h_ex = flip(h_ex);
h_ex = h_ex(1:K1+1);
H = convm(h_ex, L2).';

H = H(1:L2, 1:K1+1)

