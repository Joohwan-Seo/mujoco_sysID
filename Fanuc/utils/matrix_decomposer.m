function Mat_vec_param = matrix_decomposer(V, idx_to_delete)
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
    
    num_param = 60;
    jacob = sym(zeros(6,num_param));
    for k = 1 : 6
%         tic
        i = 7 - k;

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
%         toc


    end
    jacob = simplify(jacob);
    
%     idx_to_delete = [1,2,4,5,6,7,8,9,10,19,20];
    jacob(:,idx_to_delete) = [];
    
    Mat_vec_param = jacob;
end