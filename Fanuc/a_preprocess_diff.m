clear; close all; clc;
addpath('Dataset/')

dir = '/Dataset';

N = 4;
dt = 0.002;
%% Reading data, then store into the cell
dataset = cell(N,1);
 for i = 1 : N
    disp(i)
    name = [dir,'/square_PIH_fanuc_new_real_',num2str(i-1),'.csv'];
    T = readtable(name);
    
    dataset{i}.q = [T.q1, T.q2, T.q3, T.q4, T.q5, T.q6];
    dataset{i}.dq = [T.dq1, T.dq2, T.dq3, T.dq4, T.dq5, T.dq6];

    dataset{i}.torque = [T.t1, T.t2, T.t3, T.t4, T.t5, T.t6];
 end

acc_lim = 2;

%% differentiate to make acceleration signal
% #NOTE If the signal is noisy, we may need to use kinematic Kalman Filtering
% to first cut off the noise and to obtain the acceleration signal

for i = 1 : N
%     acc = zeros(M,2);

    acc1 = diff(dataset{i}.dq(:,1))/dt;
    acc2 = diff(dataset{i}.dq(:,2))/dt;
    acc3 = diff(dataset{i}.dq(:,3))/dt;
    acc4 = diff(dataset{i}.dq(:,4))/dt;
    acc5 = diff(dataset{i}.dq(:,5))/dt;
    acc6 = diff(dataset{i}.dq(:,6))/dt;

    acc = [acc1, acc2, acc3, acc4, acc5, acc6];


    dataset{i}.ddq = acc(500:end,:);
    dataset{i}.dq = dataset{i}.dq(500:end-1,:);
    dataset{i}.q = dataset{i}.q(500:end-1,:);
    dataset{i}.torque = dataset{i}.torque(500:end-1,:);

    disp(i)

    if i == 1
        dataset{i}.ddq(3616:3800,:) = [];
        dataset{i}.dq(3616:3800,:) = [];
        dataset{i}.q(3616:3800,:) = [];
        dataset{i}.torque(3616:3800,:) = [];
    end
    

end

%%
save('matfiles/processed_data_square_PIH.mat','dataset');