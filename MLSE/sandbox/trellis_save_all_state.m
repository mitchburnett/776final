% Viterbi algorithm
clearvars;
% define data
N = 6;
n = 0:N-1;
M = 2;
L = 2;
Nstates = M^L;
symbols = [1; -1];%; 3; -3];

% data sent
In = [1, -1, 1, 1, -1, 1];

% correlation coefficients
x = [1/8, 1/2, 5/4, 1/2, 1/8];

y = conv(In,x);
y = y(L+1:end-L);

% redefine x for convenience.
x = [5/4, 1/2, 1/8].';

% initialize
k = 1;
cur_state = [zeros(M,L), symbols];
PBM = zeros(Nstates,N);
PBMHistory = zeros(M,N);

Ik      = cur_state(:, k+L);
Ik_prev = cur_state(:, k:(k+L-1));

BM = Ik.*(2*y(k) - Ik.*x(1) - 2*Ik_prev*x(2:end));
PBMHistory(:,k) = BM;

BM = repmat(BM, [1,M]);
PBM(:,k) = reshape(BM',[Nstates,1]);

%%
for k = 2:N
    next_state = zeros(Nstates, length(cur_state) + 1);
    for m = 1:M
        tmp = [repmat(cur_state(m,:),[M,1]), symbols];
        next_state((m-1)*M+1:m*M, :) = tmp;
    end

    BM = zeros(Nstates,1);
    for l = 1:Nstates
        state = next_state(l,:);
        Ik = state(end);
        Ik_prev = state(end-1:-1:end-L);
        BM(l) = Ik*(2*y(k) - Ik*x(1) - 2*Ik_prev*x(2:end));
    end

    PBM(:,k) = PBM(:,k-1) + BM;

    % upadte
    % The first M elements of I represents the index of PBM with the
    % largest M elements.
    [Y, I] = sort(PBM(:,k), 'descend');
    
    % force consistent tree structure top-to-bottom.
    tmpI = sort(I(1:M));
    
    %cur_state = next_state(I(1:M), :);
    cur_state = next_state(tmpI, :);
    
    % save PBM for convienence for now.
    PBMHistory(:,k) = PBM(tmpI,k);

    % save the PBM and replicate to referesh the branches that have been killed
    
    
    
%    PBM = repmat(PBM(I(1:M),:), [M,1]);
    tmp = zeros(size(PBM));
    for m = 1:M
        tmp((m-1)*M+1:m*M,:) = repmat(PBM(tmpI(m),:), [M,1]);
    end
    PBM = tmp;
%      tmp = repmat(PBM(tmpI,:),[M,1]);
%      PBM(1,:) = tmp(1,:);
%      PBM(2,:) = tmp(1,:);
%      PBM(3,:) = tmp(2,:);
%      PBM(4,:) = tmp(2,:);
end









