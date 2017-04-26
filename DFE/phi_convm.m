% This file represnets my attempt to form the phi matrix as shown for the
% DFE. This is phi, or equation (34). Right now, I am sure there are better
% ways to form some of the terms. (i.e similar to the convm method of MMSE
% or how I formed matrix in 776 using convm) However, I wanted to to this
% just to make sure and fully vet the DFE and get the right answers.

%clearvars;
% load in an observation model.
load('../Models/Forney.mat');
fn = Forney.filter.'; clear Forney;
L = length(fn)-1;

% filter length of feed forward is K1+1.
K1 = 5;

% general observation model. Forney model L1 = 0, L2 = L.
L1 = 0;
L2 = L;

% phi_first I am not positive this is the best way to form the matrix. With
% the support of the matrix starting at phi(-K1,-K1) that really has thrown
% me off. I may have been able to use the convm method right away without fliping
% by just accesing the negative of the index. If this works out, it is worth seeing
% if this can be constructed in a more tractable way.
H = convm(fn,K1+1);
Hflip = flip(H, 2); %flip the matrix around vertical center.

phi_first = zeros(K1+1, K1+1);
for k = -K1:0
    for i = -K1:0
        row_start = k+K1+1;
        row_end = row_start+L2;
        % pulls the column for the first conj term in phi.
        tmp = Hflip(row_start:row_end,-i+1);
        phi_first(K1+i+1,K1+k+1) = fn*conj(tmp);
    end
end

% same goes for phi_second. I didnt even try using convm. However, it may
% be worth seeing if there is a more tractable way.
phi_second = zeros(K1+1, K1+1);
fn_tmp = [fn, zeros(1, K1)];
for k = -K1:0
    for i = -K1:0
        % notice the +2. There is an extra +1 since the second sum in phi
        % starts from l=1 and extends to L2.
        h1 = fn_tmp(-k+2:-k+L2+1);
        h2 = fn_tmp(-i+2:-i+L2+1);
        phi_second(K1+i+1,K1+k+1) = h1*h2';
    end
end

phi = phi_first-phi_second;


% Now calculating phi as shown in the book for the forney observation
% model. This is the special case of phi in the notes on DFE from class
% when L1 = 0.

fn_tmp = [zeros(1, K1), fn, zeros(1,K1-L)]; % pad fn on both sides since the sum goes to K1 this will cause it to go beyond its support.
phi_forney = zeros(K1+1, K1+1);
for k = -K1:0
    for i = -K1:0
        % this formulation the sum reaches negative indices and so the
        % added K1+1 in tmp1 and tmp2 are to center back at 0. The K1+1 in
        % the phi_forney assignment is to have the top left of phi(-K1,-K1)
        % map to (1,1).
        tmp1 = fn_tmp(K1+1:-i+K1+1);
        tmp2 = fn_tmp(i-k+K1+1:i-k+K1+1-i);
        phi_forney(K1+i+1,K1+k+1) = conj(tmp1)*tmp2.';
    end
end
        
        
        
        
        
        
        