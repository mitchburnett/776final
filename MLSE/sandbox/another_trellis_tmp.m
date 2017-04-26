% Viterbi algorithm

% define data
N = 6;
n = 0:N-1;
In = [1, -1, 1, 1, -1, 1];
% correlation coefficients
%x = [1/2, 5/4, 1/2];
x = [1/2, 1, 5/4, 1, 1/2];
L = 2;
% create sample outputs
yn = conv(In,x);
yn = yn(L+1:end-L);
% define possible states
Ik = [-1, 1];
% process variables
lambda = zeros(2,N);
branch_metric = zeros(2,1);
tmp = zeros(2,2);
% initialize
k = 1;
lambda(1,k) = Ik(1)*(2*yn(k)-x(2)*Ik(1));
lambda(2,k) = Ik(2)*(2*yn(k)-x(2)*Ik(2));
for k = 2:N
    for m = 1:2 % 1 - top node. 2 - bottom node.
        for p = 1:2 % 1 - stay on top branch. 2 - transition to bottom branch.
        branch_metric(p) = Ik(p)*(2*yn(k)-x(2)*Ik(p) - Ik(m));
        end
        tmp(:,m) = lambda(m,k-1) + branch_metric;
    end
    lambda(:,k) = max(tmp,[],2);
end
% create decision
[L, I] = max(lambda);
I_hat = Ik(I);
 