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

Ik      = cur_state(:, k+L);
Ik_prev = cur_state(:, k:(k+L-1));

PBM = Ik.*(2*y(k) - Ik.*x(1) - 2*Ik_prev*x(2:end));

k = 2;
next_state = zeros(Nstates, L+1); % where I have been plus where one for where I am going.
for m = 1:M
    tmp = [repmat(cur_state(m,2:end),[M,1]), symbols];
    next_state((m-1)*M+1:m*M, :) = tmp;
end

BM = zeros(Nstates,1);
for l = 1:Nstates
    state = next_state(l,:);
    Ik = state(end);
    Ik_prev = state(L:-1:1);
    BM(l) = Ik*(2*y(k) - Ik*x(1) - 2*Ik_prev*x(2:end));
end

tmp = repmat(PBM, [1,M]);
tmp = reshape(tmp',[Nstates,1]);

BM_tmp = tmp + BM

[Y, I] = sort(BM_tmp, 'descend');







