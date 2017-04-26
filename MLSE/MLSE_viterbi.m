% MLSE sequence estimator

% So far this is my finished MLSE esimator implementing the viterbi
% algorihtm. The several past files that I haev get me up to this code
% here.

% In this code I keep track of the state as I go. Because I have keep track
% of the max I dont think I need to traverse back through the graph.

clearvars;
load('ForneyModel.mat');
fn = f_n; clear f_n;
% define data
N = 6;
n = 0:N-1;
M = 2;
L = 2;
Nstates = M^L;
symbols = [1; -1];%; 3; -3];

% data sent
In = [1, -1, 1, 1, -1, 1].';
%idx = randi([1,4],N,1);
%In = symbols(idx);
n = 0*randn(N,1);

% correlation coefficients
x = [1/8, 1/2, 5/4, 1/2, 1/8];
%x = [fn(1), fn(2), fn(3)].';
%x = fn.';

r = In + n;
y = conv(r,x);
y = y(L+1:end-L); % Ungerboeck
%y = y(1:end-L); % Forney

% redefine x for convenience.
x = [5/4, 1/2, 1/8].';

%% Viterbi Algorithm

% Initialize state and path metric
k = 1;
cur_state = [zeros(M,L), symbols];
PM = zeros(Nstates,N);

BM = BM_Ungerboeck(cur_state, x, y(k)); % insert BM for model here.
%BM = -BM_Forney(cur_state, x, y(k));

BM = repmat(BM, [1,M]);
PM(:,k) = reshape(BM',[Nstates,1]);

for k = 2:N
    
    % create the next possible states
    [~, ydim] = size(cur_state);
    next_state = zeros(Nstates, ydim + 1);
    for m = 1:M
        tmp = [repmat(cur_state(m,:),[M,1]), symbols];
        next_state((m-1)*M+1:m*M, :) = tmp;
    end
    
    % Compute branch metric.
    BM = BM_Ungerboeck(next_state,x,y(k)); % insert BM for model here.
    %BM = -BM_Forney(next_state,x,y(k)); % FOM and POM need to subtract branch metric.

    % Update path metric and state information.
    PM(:,k) = PM(:,k-1) + BM;
    
    
    % The first M elements of I represents the index of PBM with the largest M elements.
    [Y, I] = sort(PM(:,k), 'descend'); % descend for max, ascend for min?
    % force consistent tree structure top-to-bottom. sort the sort...
    tmpI = sort(I(1:M));
    
    % save surviving paths.
    cur_state = next_state(tmpI, :);
    
    % there is no more to process. Leave final kill to decision logic.
    if k~=N
        % Determine surviving paths.
        % this way also seems to work to reform PBM
        maxBM = PM(tmpI,1:k);
        tmp = repmat(maxBM,[1,M]);
        tmp = reshape(tmp', [k,Nstates])'; %this is voodoo. Takes a while messisng with repmat/reshape to see and get this.
        PM(:,1:k) = tmp;
    end
end

% Finish the estimation and selection code. i.e. max of the last BM









