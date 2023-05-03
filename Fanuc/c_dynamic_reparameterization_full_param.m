clear;
%%
addpath('utils/')
%% Dynamic equation derivation
load('utils/Fanuc_symbolic_DH.mat');
syms q1 q2 q3 q4 q5 q6 real
syms dq1 dq2 dq3 dq4 dq5 dq6 real
syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 real
syms tau1 tau2 tau3 tau4 tau5 tau6 real

q = [q1; q2; q3; q4; q5; q6];
dq = [dq1; dq2; dq3; dq4; dq5; dq6];
ddq = [ddq1; ddq2; ddq3; ddq4; ddq5; ddq6];
tau = [tau1; tau2; tau3; tau4; tau5; tau6];

V = M * ddq + C * dq + transpose(G);

alpha = [-pi/2, pi, -pi/2, pi/2, -pi/2, pi]; % updated 4/20/2023
d = [330, 0, 0, -420, 0, -320] * 0.001; % updated 4/14/2023
a = [50, 440, -35, 0, 0, 0] * 0.001; % updated 4/14/2023
offset = [0, -pi/2, 0, 0, 0, 0]; % updated 4/20/2023

syms xx1 yy1 zz1 real
syms xx2 yy2 zz2 real
syms xx3 yy3 zz3 real
syms xx4 yy4 zz4 real
syms xx5 yy5 zz5 real
syms xx6 yy6 zz6 real

In = [xx1, yy1, zz1; 
      xx2, yy2, zz2;
      xx3, yy3, zz3;
      xx4, yy4, zz4;
      xx5, yy5, zz5;
      xx6, yy6, zz6];

syms pc_6x pc_6y pc_6z real
syms m6 real

syms pc_5x pc_5y pc_5z real
syms m5 real

syms pc_4x pc_4y pc_4z real
syms m4 real

syms pc_3x pc_3y pc_3z real
syms m3 real

syms pc_2x pc_2y pc_2z real
syms m2 real

syms pc_1x pc_1y pc_1z real
syms m1 real
  
XX6 = xx6 + m6*(pc_6y^2 + pc_6z^2);
YY6 = yy6 + m6*(pc_6x^2 + pc_6z^2);
ZZ6 = zz6 + m6*(pc_6x^2 + pc_6y^2);

XX5 = xx5 + m5*(pc_5y^2 + pc_5z^2);
YY5 = yy5 + m5*(pc_5x^2 + pc_5z^2);
ZZ5 = zz5 + m5*(pc_5x^2 + pc_5y^2);

XX4 = xx4 + m4*(pc_4y^2 + pc_4z^2);
YY4 = yy4 + m4*(pc_4x^2 + pc_4z^2);
ZZ4 = zz4 + m4*(pc_4x^2 + pc_4y^2);

XX3 = xx3 + m3*(pc_3y^2 + pc_3z^2);
YY3 = yy3 + m3*(pc_3x^2 + pc_3z^2);
ZZ3 = zz3 + m3*(pc_3x^2 + pc_3y^2);

XX2 = xx2 + m2*(pc_2y^2 + pc_2z^2);
YY2 = yy2 + m2*(pc_2x^2 + pc_2z^2);
ZZ2 = zz2 + m2*(pc_2x^2 + pc_2y^2);

XX1 = xx1 + m1*(pc_1y^2 + pc_1z^2);
YY1 = yy1 + m1*(pc_1x^2 + pc_1z^2);
ZZ1 = zz1 + m1*(pc_1x^2 + pc_1y^2);

XY1 = m1*pc_1x*pc_1y; XY2 = m2*pc_2x*pc_2y;
XY3 = m3*pc_3x*pc_3y; XY4 = m4*pc_4x*pc_4y;
XY5 = m5*pc_5x*pc_5y; XY6 = m6*pc_6x*pc_6y;

XZ1 = m1*pc_1x*pc_1z; XZ2 = m2*pc_2x*pc_2z;
XZ3 = m3*pc_3x*pc_3z; XZ4 = m4*pc_4x*pc_4z;
XZ5 = m5*pc_5x*pc_5z; XZ6 = m6*pc_6x*pc_6z;

YZ1 = m1*pc_1y*pc_1z; YZ2 = m2*pc_2y*pc_2z;
YZ3 = m3*pc_3y*pc_3z; YZ4 = m4*pc_4y*pc_4z;
YZ5 = m5*pc_5y*pc_5z; YZ6 = m6*pc_6y*pc_6z;

XX = [XX1; XX2; XX3; XX4; XX5; XX6];
YY = [YY1; YY2; YY3; YY4; YY5; YY6];
ZZ = [ZZ1; ZZ2; ZZ3; ZZ4; ZZ5; ZZ6];
XY = [XY1; XY2; XY3; XY4; XY5; XY6];
XZ = [XZ1; XZ2; XZ3; XZ4; XZ5; XZ6];
YZ = [YZ1; YZ2; YZ3; YZ4; YZ5; YZ6];

mX = [m1*pc_1x; m2*pc_2x; m3*pc_3x; m4*pc_4x; m5*pc_5x; m6*pc_6x];
mY = [m1*pc_1y; m2*pc_2y; m3*pc_3y; m4*pc_4y; m5*pc_5y; m6*pc_6y];
mZ = [m1*pc_1z; m2*pc_2z; m3*pc_3z; m4*pc_4z; m5*pc_5z; m6*pc_6z];

pc = [pc_1x, pc_1y, pc_1z;
      pc_2x, pc_2y, pc_2z;
      pc_3x, pc_3y, pc_3z;
      pc_4x, pc_4y, pc_4z;
      pc_5x, pc_5y, pc_5z;
      pc_6x, pc_6y, pc_6z;];


m = [m1; m2; m3; m4; m5; m6];
  
num_param = 10 * length(q);
jacob = sym(zeros(length(q),num_param));
%% 
for k = 1 : length(q)
    tic
    i = length(q) + 1 - k;
    
    jacob_01 = CF(V,In(i,1),1);
    V = V - XX(i)*jacob_01;
    
    jacob_02 = CF(V,In(i,2),1);
    V = V - YY(i)*jacob_02;
    
    jacob_03 = CF(V,In(i,3),1);
    V = V - ZZ(i)*jacob_03;
      
    jacob_04 = CF(V,[m(i),pc(i,1),pc(i,2)],4); %%mxy
    V = V - XY(i)*jacob_04;  
    jacob_05 = CF(V,[m(i),pc(i,1),pc(i,3)],4); %%mxz
    V = V - XZ(i)*jacob_05;  
    jacob_06 = CF(V,[m(i),pc(i,2),pc(i,3)],4); %%myz
    V = V - YZ(i)*jacob_06;
    
    jacob_07 = CF(V,[m(i),pc(i,1)],2);
    V = V - mX(i)*jacob_07;
    jacob_08 = CF(V,[m(i),pc(i,2)],2);
    V = V - mY(i)*jacob_08;
    jacob_09 = CF(V,[m(i),pc(i,3)],2);
    V = V - mZ(i)*jacob_09;
    
    V = simplify(expand(V));
    
    jacob_10 = CF(V,m(i),1);
    V = V - m(i)*jacob_10;
    
    V = simplify(expand(V));
    
    jacob(:,(i-1)*10+1:(i)*10) = ...
        [jacob_01, jacob_02, jacob_03, jacob_04, jacob_05, jacob_06,...
         jacob_07, jacob_08, jacob_09, jacob_10];
    toc
    
    
end
%%
tic
jacob = simplify(expand(jacob));
idx = find(sum(jacob,1) == 0);
toc

jacob(:,idx) = []; % removing the unnecessary indices

save('utils/idx_to_remove.mat','idx');

m = length(jacob);

fprintf('number of parameters are %d \n',m);
%% Hard coding part; adding damping and coulomb friction
jacob_53 = [dq1; 0; 0; 0; 0; 0;];
jacob_54 = [0; dq2; 0; 0; 0; 0;];
jacob_55 = [0; 0; dq3; 0; 0; 0;];
jacob_56 = [0; 0; 0; dq4; 0; 0;];
jacob_57 = [0; 0; 0; 0; dq5; 0;];
jacob_58 = [0; 0; 0; 0; 0; dq6;];

jacob_59 = -[sign(dq1); 0; 0; 0; 0; 0];
jacob_60 = -[0; sign(dq2); 0; 0; 0; 0];
jacob_61 = -[0; 0; sign(dq3); 0; 0; 0];
jacob_62 = -[0; 0; 0; sign(dq4); 0; 0];
jacob_63 = -[0; 0; 0; 0; sign(dq5); 0];
jacob_64 = -[0; 0; 0; 0; 0; sign(dq6)];

jacob_65 = [ddq1; 0; 0; 0; 0; 0;];
jacob_66 = [0; ddq2; 0; 0; 0; 0;];
jacob_67 = [0; 0; ddq3; 0; 0; 0;];
jacob_68 = [0; 0; 0; ddq4; 0; 0;];
jacob_69 = [0; 0; 0; 0; ddq5; 0;];
jacob_70 = [0; 0; 0; 0; 0; ddq6;];


% jacob = [jacob, jacob_53, jacob_54, jacob_55, jacob_56, jacob_57, jacob_58, ...
%         jacob_59, jacob_60, jacob_61, jacob_62, jacob_63, jacob_64, ...
%         jacob_65, jacob_66, jacob_67, jacob_68, jacob_69, jacob_70];

jacob = [jacob, jacob_53, jacob_54, jacob_55, jacob_56, jacob_57, jacob_58, ...
        jacob_59, jacob_60, jacob_61, jacob_62, jacob_63, jacob_64];


tic
matlabFunction(jacob,'File','utils/funJacobian_fanuc_full_param', ...
    'Vars',[q1, q2, q3, q4, q5, q6, dq1, dq2, dq3, dq4, dq5, dq6,...
            ddq1, ddq2, ddq3, ddq4, ddq5, ddq6,...
            ],'Outputs',{'jacob'});
toc
