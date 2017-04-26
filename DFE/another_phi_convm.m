% This seems to be an easier, simpler way to compute the phi matrix.
% The idea is tho think of the phi matrix as a temporal correlation matrix

% load in an observation model.
%clearvars;
load('../Models/Forney.mat');
fn = Forney.filter; clear Forney;
L = length(fn)-1;

% filter length of feed forward is K1+1.
K1 = 5;

% general observation model. Forney model L1 = 0, L2 = L.
L1 = 0;
L2 = L;

H = convm(fn,K1+1);

phi1 = zeros(K1+1, K1+1);
for k = 1:K1+1
    for i = 1:K1+1
        phi1(i,k) = H(:,i)'*H(:,k);
    end
end

H2 = convm(fn(2:end), K1+1);
phi2 = zeros(K1+1, K1+1);
for k = 1:K1+1
    for i = 1:K1+1
        phi2(i,k) = H2(L2+1:end,i)'*H2(L2+1:end,k);
    end
end

phi_tmp = phi1-phi2;
