clear; close all; clc
%%
addpath('utils/');
addpath('generated_functions/')
addpath('matfiles/')
%%
load('Fanuc_symbolic_DH.mat');
load('fanuc_result_arm.mat')

syms q1 q2 q3 q4 q5 q6 real
syms dq1 dq2 dq3 dq4 dq5 dq6 real
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 real
syms tau1 tau2 tau3 tau4 tau5 tau6 real

load('utils/idx_to_remove.mat');


transf_mat = transf(1:length(THETA)-12,1:size(transf,2)-12);
Theta = THETA(1:length(THETA)-12,1);
% Negative 12 to remove friction and damping term
%% Mass
n_state = 6;
M_reconstruct = sym(zeros(6,6));

for i = 1:n_state
    tic
    param_mat = matrix_decomposer(M(:,i), idx);
%     M_reconstruct(:,i) = param_mat * transpose(transf_mat) * Theta;
    M_reconstruct(:,i) = param_mat * transpose(transf()) * THETA;
    fprintf('Mass %d -th, ',i);
    toc    
end
M_reconstruct = M_reconstruct + diag([5, 5/2, 5/3, 5/4, 5/5, 5/6]); %% Adding Armature values

tic
matlabFunction(M_reconstruct,'File','generated_functions/funMass_UR5e', ...
    'Vars',[q1, q2, q3, q4, q5, q6,...
            ],'Outputs',{'Mass'});
toc
%% Coriolis

C_reconstruct = sym(zeros(6,6));

for i = 1:n_state
    tic
    param_mat = matrix_decomposer(C(:,i), idx);
    C_reconstruct(:,i) = param_mat * transpose(transf_mat) * Theta;
    fprintf('Coriolis %d -th, ',i);
    toc    
end

tic
matlabFunction(C_reconstruct,'File','generated_functions/funCori_UR5e', ...
    'Vars',[q1, q2, q3, q4, q5, q6, dq1, dq2, dq3, dq4, dq5, dq6...
            ],'Outputs',{'Coriolis'});
toc
%% Gravity
G_reconstruct = sym(zeros(6,1));


tic
param_mat = matrix_decomposer(transpose(G), idx);
G_reconstruct(:,1) = param_mat * transpose(transf_mat) * Theta;
fprintf('Gravity %d -th, ',1);
toc    


tic
matlabFunction(G_reconstruct,'File','generated_functions/funGrav_UR5e', ...
    'Vars',[q1, q2, q3, q4, q5, q6,...
            ],'Outputs',{'Gravity'});
toc