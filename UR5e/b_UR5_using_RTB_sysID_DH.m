clear; close all; clc;
%%
addpath('utils/')
addpath('generated_functions/')

alpha = [-pi/2, 0, 0, -pi/2, pi/2, 0];
d = [0.163, 0, 0, 0.134, 0.1, 0.098-0.004]; % d is given in z-direction
a = [0, 0.425, 0.392, 0, 0, 0]; % a is given in x-direction

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
    'alpha', alpha(1), ...
    'I', [xx1, yy1, zz1, 0, 0, 0], ... 
    'r', [pc_1x, pc_1y, pc_1z], ...       
    'm', m1); 

L(2) = Revolute('d', d(2), ...   
    'a', a(2), ...               
    'alpha', alpha(2), ...
    'I', [xx2, yy2, zz2, 0, 0, 0], ... 
    'r', [pc_2x, pc_2y, pc_2z], ...       
    'm', m2); 

L(3) = Revolute('d', d(3), ...   
    'a', a(3), ...               
    'alpha', alpha(3), ...
    'I', [xx3, yy3, zz3, 0, 0, 0], ... 
    'r', [pc_3x, pc_3y, pc_3z], ...       
    'm', m3); 

L(4) = Revolute('d', d(4), ...   
    'a', a(4), ...               
    'alpha', alpha(4), ...
    'I', [xx4, yy4, zz4, 0, 0, 0], ... 
    'r', [pc_4x, pc_4y, pc_4z], ...       
    'm', m4); 

L(5) = Revolute('d', d(5), ...   
    'a', a(5), ...               
    'alpha', alpha(5), ...
    'I', [xx5, yy5, zz5, 0, 0, 0], ... 
    'r', [pc_5x, pc_5y, pc_5z], ...       
    'm', m5); 

L(6) = Revolute('d', d(6), ...   
    'a', a(6), ...               
    'alpha', alpha(6), ...
    'I', [xx6, yy6, zz6, 0, 0, 0], ... 
    'r', [pc_6x, pc_6y, pc_6z], ...       
    'm', m6); 

% L(7) = Link('theta', pi);

UR5e = SerialLink(L, 'name', 'UR5e');

UR5e.sym()
%%
syms q1 q2 q3 q4 q5 q6
syms dq1 dq2 dq3 dq4 dq5 dq6

q = [q1,q2,q3,q4,q5,q6];

tic             
M = UR5e.inertia([q1,q2,q3,q4,q5,q6]);
C = UR5e.coriolis([q1,q2,q3,q4,q5,q6],[dq1,dq2,dq3,dq4,dq5,dq6]);
G = UR5e.gravload([q1,q2,q3,q4,q5,q6]);
toc

save('utils/ur5e_symbolic_DH.mat','M','C','G');
%%
Jb = UR5e.jacobe(q);
g_st_ = UR5e.fkine(q);
g_st = simplify([g_st_.R, g_st_.t;
                 zeros(1,3),1]);
             
Je = UR5e.jacob0(q);

Jb = simplify(Jb);
g_st = simplify(g_st);
Je = simplify(Je);
tic
matlabFunction(Jb,'File','generated_functions/Jb','Vars',[q1, q2, q3, q4, q5, q6],'Outputs',{'Js_mat'});
matlabFunction(Je,'File','generated_functions/Je','Vars',[q1, q2, q3, q4, q5, q6],'Outputs',{'Js_mat'});
matlabFunction(g_st,'File','generated_functions/g_st','Vars',[q1, q2, q3, q4, q5, q6],'Outputs',{'g_st_mat'});
toc