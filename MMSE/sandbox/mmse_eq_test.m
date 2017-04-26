% MMSE Equalizer
%clearvars;
% load in an observation model.
load('ForneyModel.mat');
fn_norm = f_n; clear f_n;
load('F_noNorm.mat');
fn = f_n; clear f_n;
load('xn_short.mat');
L = length(fn_norm)-1;

fn_tmp = [zeros(1,L), fn_norm, zeros(1,L)];
Z = L+1;
tmpx = zeros(1,2*L+1);
for k = -L:L
    for n = 0:L
        tmpx(k+Z) = tmpx(k+Z) + conj(fn_tmp(Z+n))*fn_tmp(n+k + Z);
    end
end

%%
fn_tmp2 = convm(fn_norm,2*L+1);
tmpx2 = zeros(2*L+1, 2*L+1);
for l = 1:2*L+1
    for i = 1:2*L+1
        tmpx2(l,i) = fn_tmp2(:,i)'*fn_tmp2(:,l);
    end
end

% V2 = convm(v2,p);
% r_vx_hat = zeros(p,1);
% for k = 1:p
%    r_vx_hat(k,1) = x'*V2(1:samples,k);
% end
% r_vx_hat = 1/samples*r_vx_hat;
% 
% r_v_hat = zeros(p,p);
% for k = 1:p
%     for i = 1:p
%         r_v_hat(k,i) = V2(1:samples,i)'*V2(1:samples,k);
%     end
% end