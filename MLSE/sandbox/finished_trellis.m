% MLSE sequence estimator

% So far this is my finished MLSE esimator implementing the viterbi
% algorihtm. The several past files that I haev get me up to this code
% here.

% In this code I keep track of the state as I go. Because I have keep track
% of the max I dont think I need to traverse back through the graph.

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
n = 0*randn(1,N);

% correlation coefficients
x = [1/8, 1/2, 5/4, 1/2, 1/8];

r = In + n;
y = conv(r,x);
y = y(L+1:end-L);

% redefine x for convenience.
x = [5/4, 1/2, 1/8].';

% initialize state and partial path metric
k = 1;
cur_state = [zeros(M,L), symbols];
PPM = zeros(Nstates,N);

Ik      = cur_state(:, k+L);
Ik_prev = cur_state(:, k:(k+L-1));

BM = Ik.*(2*y(k) - Ik.*x(1) - 2*Ik_prev*x(2:end));


BM = repmat(BM, [1,M]);
PPM(:,k) = reshape(BM',[Nstates,1]);

%%
for k = 2:N
    
    % create the next possible states
    [~, ydim] = size(cur_state);
    next_state = zeros(Nstates, ydim + 1);
    for m = 1:M
        tmp = [repmat(cur_state(m,:),[M,1]), symbols];
        next_state((m-1)*M+1:m*M, :) = tmp;
    end
    
    % Compute branch metric.
%     BM = zeros(Nstates,1);
%     for l = 1:Nstates
%         state = next_state(l,:);
%         BM(l) = BM_Ungerboeck(state, x, y(k));
%     end
    BM = BM_Ungerboeck(next_state,x,y(k));

    % Update partial branch metric and state information.
    PPM(:,k) = PPM(:,k-1) + BM;
    
    % The first M elements of I represents the index of PBM with the largest M elements.
    [Y, I] = sort(PPM(:,k), 'descend');
    
    % force consistent tree structure top-to-bottom.
    tmpI = sort(I(1:M));
    
    cur_state = next_state(tmpI, :);
    
    % there is no more to process. Leave final kill to decision logic.
    if k~=N
        % Determine surviving paths.
        % this way also seems to work to reform PBM
        maxBM = PPM(tmpI,1:k);
        tmp = repmat(maxBM,[1,M]);
        tmp = reshape(tmp', [k,Nstates])'; %this is voodoo. Takes a while messisng with repmat/reshape to see and get this.
        PPM(:,1:k) = tmp;
    end
end









