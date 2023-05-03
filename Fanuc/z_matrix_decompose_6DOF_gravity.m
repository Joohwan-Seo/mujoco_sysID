clear;
%%
addpath('utils/');
addpath('generated_functions/')
addpath('matfiles/')
%%
load('Fanuc_symbolic_DH.mat');
load('fanuc_result_arm_gravity.mat')

syms q1 q2 q3 q4 q5 q6 real
syms dq1 dq2 dq3 dq4 dq5 dq6 real
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 real
syms tau1 tau2 tau3 tau4 tau5 tau6 real

load('utils/idx_to_remove_grav.mat');


transf_mat = transf(1:length(THETA),1:size(transf,2));
Theta = THETA(1:length(THETA),1);
% Negative 12 to remove friction and damping term
%% Gravity
G_reconstruct = sym(zeros(6,1));


tic
param_mat = matrix_decomposer(transpose(G), idx);
G_reconstruct(:,1) = param_mat * transpose(transf_mat) * Theta;
fprintf('Gravity %d -th, ',1);
toc    


tic
matlabFunction(G_reconstruct,'File','generated_functions/funGrav_Fanuc', ...
    'Vars',[q1, q2, q3, q4, q5, q6,...
            ],'Outputs',{'Gravity'});
toc