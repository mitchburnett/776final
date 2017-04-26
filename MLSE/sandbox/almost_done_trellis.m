% Viterbi algorithm
clearvars;
% define data
N = 6;
n = 0:N-1;
M = 4;
L = 2;
Nstates = M^L;
symbols = [1; -1; 3; -3];

% data sent
In = [1, -3, 1, 3, -1, 1];
n = .5*randn(1,N);

% correlation coefficients
x = [1/8, 1/2, 5/4, 1/2, 1/8];

r = In + n;
y = conv(r,x);
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
    [~, ydim] = size(cur_state);
    next_state = zeros(Nstates, ydim + 1);
    for m = 1:M
        tmp = [repmat(cur_state(m,:),[M,1]), symbols];
        next_state((m-1)*M+1:m*M, :) = tmp;
    end
    
    % Compute branch metric.
    BM = zeros(Nstates,1);
    for l = 1:Nstates
        state = next_state(l,:);
        Ik = state(end);
        Ik_prev = state(end-1:-1:end-L);
        BM(l) = Ik*(2*y(k) - Ik*x(1) - 2*Ik_prev*x(2:end));
    end

    % Update partial branch metric and state information.
    PBM(:,k) = PBM(:,k-1) + BM;
    
    % The first M elements of I represents the index of PBM with the largest M elements.
    [Y, I] = sort(PBM(:,k), 'descend');
    
    % force consistent tree structure top-to-bottom.
    tmpI = sort(I(1:M));
    
    %cur_state = next_state(I(1:M), :);
    cur_state = next_state(tmpI, :);
    
    % save PBM for convienence for now.
    PBMHistory(:,k) = PBM(tmpI,k);
    
    % there is no more to process. Leave final kill to decision logic.
    if k~=N
        % this way also seems to work to reform PBM
        maxBM = PBM(tmpI,1:k);
        tmp = repmat(maxBM,[1,M]);
        tmp = reshape(tmp', [k,Nstates])'; %this is voodoo. Takes a while messisng with repmat/reshape to see and get this.
        PBM(:,1:k) = tmp;
    end
end









