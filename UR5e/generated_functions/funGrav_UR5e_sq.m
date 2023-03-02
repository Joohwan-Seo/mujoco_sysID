function Gravity = funGrav_UR5e(q1,q2,q3,q4,q5,q6)
%FUNGRAV_UR5E
%    GRAVITY = FUNGRAV_UR5E(Q1,Q2,Q3,Q4,Q5,Q6)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    16-Jan-2023 21:41:41

t2 = cos(q2);
t3 = cos(q3);
t4 = cos(q4);
t5 = cos(q5);
t6 = cos(q6);
t7 = sin(q2);
t8 = sin(q3);
t9 = sin(q4);
t10 = sin(q5);
t11 = sin(q6);
t12 = q2+q3;
t16 = -q5;
t13 = cos(t12);
t14 = q4+t12;
t15 = sin(t12);
t17 = cos(t14);
t18 = q5+t14;
t19 = sin(t14);
t22 = t14+t16;
t25 = t15.*3.877312284323098e-5;
t28 = t6.*t9.*t13.*1.141457811669069e-4;
t29 = t4.*t6.*t15.*1.141457811669069e-4;
t30 = t9.*t11.*t13.*4.431809231369314e-5;
t31 = t4.*t11.*t15.*4.431809231369314e-5;
t36 = t4.*t5.*t11.*t13.*1.141457811669069e-4;
t37 = t4.*t5.*t6.*t13.*4.431809231369314e-5;
t38 = t5.*t9.*t11.*t15.*1.141457811669069e-4;
t39 = t5.*t6.*t9.*t15.*4.431809231369314e-5;
t46 = t13.*1.639881848811534e+1;
t20 = cos(t18);
t21 = sin(t18);
t23 = cos(t22);
t24 = sin(t22);
t26 = -t25;
t27 = t17.*6.308654300135727e-5;
t33 = -t28;
t34 = -t29;
t40 = -t36;
t41 = -t37;
t44 = t10.*t17.*2.505059027848395e-1;
t45 = t19.*1.871665046617155;
t47 = -t46;
t32 = t20.*4.680505562566946e-6;
t42 = t23.*4.680505562566946e-6;
t48 = t21.*3.782291805132319e-1;
t50 = t24.*3.782291805132319e-1;
t35 = -t32;
t43 = -t42;
t49 = -t48;
et1 = t2.*(-4.001833845165367e+1)-t7.*2.306711275362242e-4;
et2 = t21.*(-2.529762291208122e-1)+t24.*2.529762291208122e-1+t26+t27+t35+t43+t45+t47;
et3 = t2.*t3.*t6.*t9.*(-1.141457811669069e-4)-t2.*t4.*t6.*t8.*1.141457811669069e-4-t3.*t4.*t6.*t7.*1.141457811669069e-4+t2.*t3.*t9.*t11.*4.431809231369314e-5;
et4 = t2.*t4.*t8.*t11.*4.431809231369314e-5+t3.*t4.*t7.*t11.*4.431809231369314e-5+t6.*t7.*t8.*t9.*1.141457811669069e-4-t7.*t8.*t9.*t11.*4.431809231369314e-5;
et5 = t2.*t3.*t4.*t5.*t6.*(-4.431809231369314e-5)-t2.*t3.*t4.*t5.*t11.*1.141457811669069e-4+t2.*t5.*t6.*t8.*t9.*4.431809231369314e-5+t3.*t5.*t6.*t7.*t9.*4.431809231369314e-5;
et6 = t4.*t5.*t6.*t7.*t8.*4.431809231369314e-5+t2.*t5.*t8.*t9.*t11.*1.141457811669069e-4+t3.*t5.*t7.*t9.*t11.*1.141457811669069e-4+t4.*t5.*t7.*t8.*t11.*1.141457811669069e-4;
mt1 = [0.0,et1+et2+et3+et4+et5+et6,t26+t27+t30+t31+t33+t34+t35+t38+t39+t40+t41+t43+t44+t45+t47+t49+t50,t27+t30+t31+t33+t34+t35+t38+t39+t40+t41+t43+t44+t45+t49+t50];
mt2 = [t35+t42+t49-t50+t5.*t19.*2.505059027848395e-1+t6.*t10.*t19.*4.431809231369314e-5+t10.*t11.*t19.*1.141457811669069e-4];
mt3 = [t6.*t17.*(-4.431809231369314e-5)-t11.*t17.*1.141457811669069e-4-t5.*t6.*t19.*1.141457811669069e-4+t5.*t11.*t19.*4.431809231369314e-5];
Gravity = reshape([mt1,mt2,mt3],6,1);