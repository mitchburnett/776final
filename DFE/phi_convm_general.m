% This seems to be an easier, simpler way to compute the phi matrix.
% The idea is tho think of the phi matrix as a temporal correlation matrix
%% Load in an observation model -- Options ['Ungerboeck', 'Forney', 'PSMFOM']
model = 'PSMFOM';

load(['../Models/', model, '.mat']);
L1 = eval([model, '.L1']);
L2 = eval([model, '.L2']);
hn = eval([model, '.filter']);
eval(['clear ' model]);

% filter length of feed forward is K1+1.
K1 = 5;
K2 = L2;
% % general observation model. Forney model L1 = 0, L2 = L.
% L1 = 0;
% L2 = L;

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
