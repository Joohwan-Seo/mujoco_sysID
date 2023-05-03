clear;
%%
addpath('matfiles/')
addpath('utils/')
addpath('generated_functions/')

load('matfiles/processed_data_square_PIH.mat');

force_stack = [];
W_stack = [];

N = 4;
init = 1;
dummy = cell(N,1);
%%
Tau_set = dataset{1}.torque;
Tau_set = [Tau_set; dataset{2}.torque];
Tau_set = [Tau_set; dataset{3}.torque];
Tau_set = [Tau_set; dataset{4}.torque];
%%
% std_tau = [5.8901   40.5046   13.0007    1.1680    0.2345    0.0559];
% std_tau = [5.8901   40.5046   13.0007    1.1680    0.2345    0.001];
% std_tau = [1.8061,1.3445,1.4437,1.1833,2.1091,1.8721];
% std_tau = std(Tau_set);

% std_tau = [1.7341   20.8867   26.7979   13.2446    6.5856    1.4974];

std_tau = ones(1,6);

grav = @funGrav_Fanuc;
%%
for i = init : N
    disp(i)
    q = dataset{i}.q;
    dq = dataset{i}.dq;
    ddq = dataset{i}.ddq;
    tau = dataset{i}.torque;
    if iscell(tau)
    end
    tic
    M = length(q);
    force_stack = zeros(6*M,1);
%     force_stack_collection = zeros(M,6);
%     W_stack = zeros(6*M,70);
    W_stack = zeros(6*M,64);
    for k = 1 : M
        armature_term = [5*ddq(k,1), 5/2*ddq(k,2), 5/3*ddq(k,3), 5/4*ddq(k,4), 5/5*ddq(k,5), 5/6*ddq(k,6)];
        force_stack(6*(k-1)+1:6*k) = (((tau(k,:) - armature_term))' - grav(q(k,1),q(k,2),q(k,3),q(k,4),q(k,5),q(k,6)))./std_tau';
%         force_stack_collection(k,:) = tau(k,:) - armature_term;
%         force_stack(6*(k-1)+1:6*k) = ((tau(k,:))./std_tau)';
%         jacob_part = funJacobian_full_param_fanuc(q(k,1),q(k,2),q(k,3),q(k,4),q(k,5),q(k,6),...
%                                   dq(k,1),dq(k,2),dq(k,3),dq(k,4),dq(k,5),dq(k,6),...
%                                   ddq(k,1),ddq(k,2),ddq(k,3),ddq(k,4),ddq(k,5),ddq(k,6));
        jacob_part = funJacobian_fanuc_wo_grav(q(k,1),q(k,2),q(k,3),q(k,4),q(k,5),q(k,6),...
                                  dq(k,1),dq(k,2),dq(k,3),dq(k,4),dq(k,5),dq(k,6),...
                                  ddq(k,1),ddq(k,2),ddq(k,3),ddq(k,4),ddq(k,5),ddq(k,6));
        W_stack(6*(k-1)+1:6*k, :) = jacob_part./std_tau';
    end
    dummy{i}.W_stack = W_stack;
    dummy{i}.force_stack = force_stack;
%     dummy{i}.force_stack_collection = force_stack_collection;
    toc
end

W_big_stack = [];
force_big_stack = [];
force_collection_stack = [];
for i = init : N
    W_big_stack = [W_big_stack; dummy{i}.W_stack];
    force_big_stack = [force_big_stack; dummy{i}.force_stack];
%     force_collection_stack = [force_stack_collection; dummy{i}.force_stack_collection];
end

rank_W = rank(W_big_stack'*W_big_stack);
fprintf('Rank of W^T*W matrix is :%d.\n',rank_W);
%%
W_full = rref(W_big_stack);
transf = W_full(1:rank_W,:);

idxs = find(abs(transf)<1e-08);
transf(idxs) = 0;

W_big_stack_reduced = W_big_stack * transf'; 

rank_reduced = rank(W_big_stack_reduced'*W_big_stack_reduced);

fprintf('Rank of reduced regressor matrix is :%d.\n',rank_reduced);
fprintf('Confirm that two values above are identical!')
%%
% idx = find(sig_W_big == 0);

sig_W_big = std(W_big_stack_reduced,1);

W_big_stack_norm = (W_big_stack_reduced)./sig_W_big;

% sig_W_big(:,idx) = [];
% mean_W_big(:,idx) = [];
%%
mdl = fitlm(W_big_stack_reduced,force_big_stack, 'intercept',false);
THETA = mdl.Coefficients.Estimate;

%% final cost
J = norm(force_big_stack - W_big_stack_reduced * THETA,2)

force_recovered = W_big_stack_reduced * THETA;
%%
N = length(force_big_stack)/6;
tau = reshape(force_big_stack,[6,N]);

tau_r = reshape(force_recovered,[6,N]);

figure
plot(tau(1,:)); hold on; plot(tau_r(1,:));
legend('original force','recovered force')
title('Link 1')

figure
plot(tau(2,:)); hold on; plot(tau_r(2,:));
legend('original force','recovered force')
title('Link 2')

figure
plot(tau(3,:)); hold on; plot(tau_r(3,:));
legend('original force','recovered force')
title('Link 3')

figure
plot(tau(4,:)); hold on; plot(tau_r(4,:));
legend('original force','recovered force')
title('Link 4')

figure
plot(tau(5,:)); hold on; plot(tau_r(5,:));
legend('original force','recovered force')
title('Link 5')

figure
plot(tau(6,:)); hold on; plot(tau_r(6,:));
legend('original force','recovered force')
title('Link 6')

%%
save('matfiles/fanuc_result_wo_grav.mat','THETA', 'transf', 'rank_W');