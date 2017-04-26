% need to save the most probably for each ending state not just the M most
% probable.

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
M = 4;
L = 3;
Nstates = M^L;
%symbols = [-1-1j; -1+1j; 1-1j; 1+1j];
symbols = [1; -1; 3; -3];

% data sent
In = [1, -3, 1, 3, -1, 1].';
%idx = randi([1,4],N,1);
%In = symbols(idx);
n = 0*randn(N,1);

% correlation coefficients
x = [1/32, 1/8, 1/2, 5/4, 1/2, 1/8, 1/32];
%x = [fn(1), fn(2), fn(3)].';
%x = fn.';

r = In + n;
y = conv(r,x);
y = y(L+1:end-L); % Ungerboeck
%y = y(1:end-L); % Forney

% redefine x for convenience.
x = [5/4, 1/2, 1/8, 1/32].';

%% Viterbi Algorithm

% Initialize state and path metric
k = 1; 
%cur_state = [zeros(M^L,L), repmat(symbols, [M^(L-1),1])];
next_state = [zeros(M^L,L), repmat(symbols, [M^(L-1),1])];
cur_state = next_state;
PM = zeros(Nstates,N);

BM = BM_Ungerboeck(next_state, x, y(k));
%BM = -BM_Forney(cur_state, x, y(k));

PM(:,k) = BM;

tmpPM = [];
    for m = 1:M^(L-1)
        % Group all the states for a specific state ending
        iso_state = next_state(m:M^(L-1)/M:end,:);
        iso_PM = PM(m:M^(L-1)/M:end,:);
        [X, idx] = sort(iso_PM(:,k));
        tmpIdx = sort(idx(1:M)); % this might be problematic in the future (I mean to say I am not sure if M is the dynamic way to choose how many to save.)
        
        cur_state((m-1)*M+1:m*M,:) = iso_state(tmpIdx, :);
        
        maxBM = iso_PM(tmpIdx, 1:k);
        
        % update the PM. store as a tmp as to not overwrite until finished through original PM.
        tmpPM = [tmpPM; maxBM];
    end
    
PM(:,k) = tmpPM;

for k = 2:N
    
    % create the next possible states
    [~, ydim] = size(cur_state);
    next_state = zeros(Nstates, ydim + 1);
    
    next_state(:,1:ydim) = cur_state;
    next_state(:,end) = repmat(symbols, [M^(L-1),1]);
    cur_state = next_state;
    
    % Compute branch metric.
    BM = BM_Ungerboeck(next_state,x,y(k)); % insert BM for model here.
    %BM = -BM_Forney(next_state,x,y(k)); % FOM and POM need to subtract branch metric.

    % Update path metric and state information.
    PM(:,k) = PM(:,k-1) + BM;
    
    %%%
    % added code to save most probable branch of each possible state
    % ending.
    tmpPM = [];
    for m = 1:M
        % Group all the states for a specific state ending
        iso_state = next_state(m:M^(L-1)/M:end,:);
        iso_PM = PM(m:M^(L-1)/M:end,:);
        [X, idx] = sort(iso_PM(:,k), 'descend');
        tmpIdx = sort(idx(1:M^(L-1))); % this might be problematic in the future (I mean to say I am not sure if M is the dynamic way to choose how many to save.)
        
        % If I build next state here I need to build the first next state
        % outside the loop. i.e make cur_state the same size as next_state.
        % I also need to be careful as to build the cur_state correctly as
        % to keep the same sorted order to line up PM additions.
        cur_state((m-1)*M^(L-1)+1:m*M^(L-1),:) = iso_state(tmpIdx, :);
        
        maxBM = iso_PM(tmpIdx, k);
        
        % update the PM. store as a tmp as to not overwrite until finished
        % through original PM.
        tmpPM = [tmpPM; maxBM];
    end
    
    PM(:,k) = tmpPM;
    
end

% Finish the estimation and selection code. i.e. max of the last BM









