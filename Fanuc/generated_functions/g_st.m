function g_st_mat = g_st(q1,q2,q3,q4,q5,q6)
%G_ST
%    G_ST_MAT = G_ST(Q1,Q2,Q3,Q4,Q5,Q6)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    21-Apr-2023 12:25:53

t2 = cos(q1);
t3 = cos(q2);
t4 = cos(q3);
t5 = cos(q4);
t6 = cos(q5);
t7 = cos(q6);
t8 = sin(q1);
t9 = sin(q2);
t10 = sin(q3);
t11 = sin(q4);
t12 = sin(q5);
t13 = sin(q6);
t14 = -q3;
t17 = pi./2.0;
t15 = t2.*t5;
t16 = t5.*t8;
t18 = q2+t14;
t19 = -t17;
t22 = t2.*t6.*t11;
t23 = t6.*t8.*t11;
t25 = t2.*t3.*t10.*t11;
t26 = t2.*t4.*t9.*t11;
t27 = t3.*t8.*t10.*t11;
t28 = t4.*t8.*t9.*t11;
t20 = cos(t18);
t21 = sin(t18);
t24 = q2+t19;
t30 = -t22;
t33 = t3.*t6.*t10.*t15;
t34 = t4.*t6.*t9.*t15;
t35 = t3.*t6.*t10.*t16;
t36 = t4.*t6.*t9.*t16;
t37 = -t25;
t38 = -t28;
t29 = cos(t24);
t31 = sin(t24);
t32 = t12.*t21;
t39 = t5.*t6.*t20;
t40 = t2.*t12.*t20;
t41 = t8.*t12.*t20;
t42 = -t34;
t43 = -t36;
t45 = t16+t26+t37;
t46 = t15+t27+t38;
t44 = t32+t39;
t47 = t23+t33+t40+t42;
t48 = t30+t35+t41+t43;
mt1 = [-t7.*t47-t13.*t45,t13.*t46+t7.*(t22-t35+t36-t41),t7.*t44-t11.*t13.*t20,0.0,t7.*t45-t13.*t47,-t7.*t46+t13.*(t22-t35+t36-t41),t13.*t44+t7.*t11.*t20,0.0,-t12.*(t8.*t11-t5.*(t2.*t4.*t29+t2.*t10.*t31))-t6.*(t2.*t4.*t31-t2.*t10.*t29),t12.*(t2.*t11-t3.*t10.*t16+t4.*t9.*t16)+t6.*t8.*t20,-t6.*t21+t5.*t12.*t20,0.0];
mt2 = [t2./2.0e+1+t2.*t9.*(1.1e+1./2.5e+1)+t2.*t3.*t4.*(2.1e+1./5.0e+1)-t2.*t3.*t10.*(7.0./2.0e+2)+t2.*t4.*t9.*(7.0./2.0e+2)+t2.*t9.*t10.*(2.1e+1./5.0e+1)-t8.*t11.*t12.*(8.0./2.5e+1)+t2.*t3.*t4.*t6.*(8.0./2.5e+1)+t2.*t6.*t9.*t10.*(8.0./2.5e+1)-t3.*t10.*t12.*t15.*(8.0./2.5e+1)+t4.*t9.*t12.*t15.*(8.0./2.5e+1)];
mt3 = [t8./2.0e+1+t8.*t9.*(1.1e+1./2.5e+1)+t3.*t4.*t8.*(2.1e+1./5.0e+1)-t3.*t8.*t10.*(7.0./2.0e+2)+t4.*t8.*t9.*(7.0./2.0e+2)+t2.*t11.*t12.*(8.0./2.5e+1)+t8.*t9.*t10.*(2.1e+1./5.0e+1)+t3.*t4.*t6.*t8.*(8.0./2.5e+1)+t6.*t8.*t9.*t10.*(8.0./2.5e+1)-t3.*t10.*t12.*t16.*(8.0./2.5e+1)+t4.*t9.*t12.*t16.*(8.0./2.5e+1),t3.*(1.1e+1./2.5e+1)+t3.*t4.*(7.0./2.0e+2)+t3.*t10.*(2.1e+1./5.0e+1)-t4.*t9.*(2.1e+1./5.0e+1)+t9.*t10.*(7.0./2.0e+2)+t3.*t6.*t10.*(8.0./2.5e+1)-t4.*t6.*t9.*(8.0./2.5e+1)+t3.*t4.*t5.*t12.*(8.0./2.5e+1)+t5.*t9.*t10.*t12.*(8.0./2.5e+1)+3.3e+1./1.0e+2,1.0];
g_st_mat = reshape([mt1,mt2,mt3],4,4);
