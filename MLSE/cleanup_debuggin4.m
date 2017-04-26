% MLSE sequence estimator

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

tmpPM = zeros(Nstates,1);
for m = 1:M
    % Group all the states for a specific state ending
    iso_state = next_state(m:M^(L-1)/M:end,:);
    iso_PM = PM(m:M^(L-1)/M:end,:);
    [X, idx] = sort(iso_PM(:,k), 'descend');
    tmpIdx = sort(idx(1:M^(L-1))); % this might be problematic in the future (I mean to say I am not sure if M is the dynamic way to choose how many to save.)

    cur_state((m-1)*M^(L-1)+1:m*M^(L-1),:) = iso_state(tmpIdx, :);

    maxBM = iso_PM(tmpIdx, k);

    % update the PM. store as a tmp as to not overwrite until finished through original PM.
    tmpPM((m-1)*M^(L-1)+1:m*M^(L-1)) = maxBM;
end
    
PM(:,k) = tmpPM;

%%
k = 2;

% create the next possible states
[~, ydim] = size(cur_state);
next_state = zeros(Nstates, ydim + 1);


next_state(:,1:ydim) = cur_state;
next_state(:,end) = repmat(symbols, [M^(L-1),1]);
cur_state = next_state;

BM = BM_Ungerboeck(next_state, x, y(k));
%BM = -BM_Forney(cur_state, x, y(k));

PM(:,k) = PM(:,k-1) + BM;

tmpPM = zeros(Nstates,1);
for m = 1:M
    % Group all the states for a specific state ending
    iso_state = next_state(m:M^(L-1)/M:end,:);
    iso_PM = PM(m:M^(L-1)/M:end,:);
    [X, idx] = sort(iso_PM(:,k), 'descend');
    tmpIdx = sort(idx(1:M^(L-1))); % this might be problematic in the future (I mean to say I am not sure if M is the dynamic way to choose how many to save.)

    cur_state((m-1)*M^(L-1)+1:m*M^(L-1),:) = iso_state(tmpIdx, :);

    maxBM = iso_PM(tmpIdx, k);

    % update the PM. store as a tmp as to not overwrite until finished through original PM.
    tmpPM((m-1)*M^(L-1)+1:m*M^(L-1)) = maxBM;
end

PM(:,k) = tmpPM;

%%

for k = 3:N
    
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
    
    tmpPM = zeros(Nstates,1);
    tmpState = zeros(M^(L-1), ydim + 1);
    for m = 1:M^(L-1)
        % Group all the states for a specific state ending
        iso_state = next_state((m-1)*M+1:m*M,:)
        iso_PM = PM((m-1)*M+1:m*M,:);
        [X, idx] = sort(iso_PM(:,k), 'descend');
        %tmpIdx = sort(idx); % this might be problematic in the future (I mean to say I am not sure if M is the dynamic way to choose how many to save.)
        
        %cur_state((m-1)*M^(L-1)+1:m*M^(L-1),:) = iso_state(tmpIdx, :);
        tmpState(m,:) = iso_state(idx(1),:);
        
        maxBM = iso_PM(idx(1),k);
        %maxBM = iso_PM(tmpIdx, k);
        
        % update the PM. store as a tmp as to not overwrite until finished through original PM.
        tmpPM((m-1)*M+1:m*M) = maxBM;
    end
    
    PM(:,k) = tmpPM;
    
end

% Finish the estimation and selection code. i.e. max of the last BM









