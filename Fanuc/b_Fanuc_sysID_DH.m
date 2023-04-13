clear; close all; clc;
%%
% addpath('utils/')
addpath('generated_functions/')

alpha = [-pi/2, pi, pi/2, -pi/2, pi/2, pi];
d = [330, 0, 0, -420, 0, -230] * 0.001; % d is given in z-direction
a = [50, 440, -35, 0, 0, -17.5] * 0.001; % a is given in x-direction
offset = [0, -pi/2, pi, 0, 0, pi];

syms xx1 yy1 zz1
syms xx2 yy2 zz2
syms xx3 yy3 zz3
syms xx4 yy4 zz4
syms xx5 yy5 zz5
syms xx6 yy6 zz6

syms pc_6x pc_6y pc_6z
syms m6

syms pc_5x pc_5y pc_5z
syms m5

syms pc_4x pc_4y pc_4z
syms m4

syms pc_3x pc_3y pc_3z
syms m3

syms pc_2x pc_2y pc_2z
syms m2

syms pc_1x pc_1y pc_1z
syms m1

%%

L(1) = Revolute('d', d(1), ...   
    'a', a(1), ...        
    'offset', offset(1), ...
    'alpha', alpha(1), ...
    'I', [xx1, yy1, zz1, 0, 0, 0], ... 
    'r', [pc_1x, pc_1y, pc_1z], ...       
    'm', m1); 

L(2) = Revolute('d', d(2), ...   
    'a', a(2), ...
    'offset', offset(2), ...
    'alpha', alpha(2), ...
    'I', [xx2, yy2, zz2, 0, 0, 0], ... 
    'r', [pc_2x, pc_2y, pc_2z], ...       
    'm', m2); 

L(3) = Revolute('d', d(3), ...   
    'a', a(3), ...           
    'offset', offset(3), ...
    'alpha', alpha(3), ...
    'I', [xx3, yy3, zz3, 0, 0, 0], ... 
    'r', [pc_3x, pc_3y, pc_3z], ...       
    'm', m3); 

L(4) = Revolute('d', d(4), ...   
    'a', a(4), ...      
    'offset', offset(4), ...
    'alpha', alpha(4), ...
    'I', [xx4, yy4, zz4, 0, 0, 0], ... 
    'r', [pc_4x, pc_4y, pc_4z], ...       
    'm', m4); 

L(5) = Revolute('d', d(5), ...   
    'a', a(5), ...        
    'offset', offset(5), ...
    'alpha', alpha(5), ...
    'I', [xx5, yy5, zz5, 0, 0, 0], ... 
    'r', [pc_5x, pc_5y, pc_5z], ...       
    'm', m5); 

L(6) = Revolute('d', d(6), ...   
    'a', a(6), ...   
    'offset', offset(6), ...
    'alpha', alpha(6), ...
    'I', [xx6, yy6, zz6, 0, 0, 0], ... 
    'r', [pc_6x, pc_6y, pc_6z], ...       
    'm', m6); 

% L(7) = Link('theta', pi);

Fanuc = SerialLink(L, 'name', 'Fanuc_LR_MATE_200ID');

Fanuc.sym()
%%
syms q1 q2 q3 q4 q5 q6
syms dq1 dq2 dq3 dq4 dq5 dq6

q = [q1,q2,q3,q4,q5,q6];

tic             
M = Fanuc.inertia([q1,q2,q3,q4,q5,q6]);
C = Fanuc.coriolis([q1,q2,q3,q4,q5,q6],[dq1,dq2,dq3,dq4,dq5,dq6]);
G = Fanuc.gravload([q1,q2,q3,q4,q5,q6]);
toc

save('utils/Fanuc_symbolic_DH.mat','M','C','G');
%%
Jb = Fanuc.jacobe(q);
g_st_ = Fanuc.fkine(q);
g_st = simplify([g_st_.R, g_st_.t;
                 zeros(1,3),1]);
             
Je = Fanuc.jacob0(q);

Jb = simplify(Jb);
g_st = simplify(g_st);
Je = simplify(Je);

tic
matlabFunction(Jb,'File','generated_functions/Jb','Vars',[q1, q2, q3, q4, q5, q6],'Outputs',{'Js_mat'});
matlabFunction(Je,'File','generated_functions/Je','Vars',[q1, q2, q3, q4, q5, q6],'Outputs',{'Js_mat'});
matlabFunction(g_st,'File','generated_functions/g_st','Vars',[q1, q2, q3, q4, q5, q6],'Outputs',{'g_st_mat'});
toc