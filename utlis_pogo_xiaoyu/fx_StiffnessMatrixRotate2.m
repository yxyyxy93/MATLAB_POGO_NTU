function SM_rotated = fx_StiffnessMatrixRotate2(SM_origin,A,B,C)
% A,B,C are the Euler angles under ZXZ convention  (Bunge convention), unit: rad
% A,B,C refer to the rotation of Coordinate Axes instead of the crystal
% lattice
% SM_x is the stiffness matrix

T1=[    cos(A), sin(A), 0;
        -sin(A),cos(A), 0;
        0,      0,      1];
T2=[    1,  0,      0;
        0,  cos(B), sin(B);
        0,  -sin(B),cos(B)];
T3=[    cos(C), sin(C), 0;
        -sin(C),cos(C), 0;
        0,      0,      1];

T=T3*T2*T1;
l1=T(1,1);
l2=T(2,1);
l3=T(3,1);
m1=T(1,2);
m2=T(2,2);
m3=T(3,2);
n1=T(1,3);
n2=T(2,3);
n3=T(3,3);

T_strain=[      l1^2,	m1^2,   n1^2,   2*m1*n1,          2*n1*l1,          2*l1*m1;
                l2^2,	m2^2,   n2^2,   2*m2*n2,          2*n2*l2,          2*l2*m2;
                l3^2,	m3^2,   n3^2,   2*m3*n3,          2*n3*l3,          2*l3*m3;
                l2*l3, m2*m3, n2*n3,m2*n3+m3*n2,    n2*l3+n3*l2,    l2*m3+l3*m2;
                l1*l3, m1*m3, n1*n3,m1*n3+m3*n1,    n1*l3+n3*l1,    l1*m3+l3*m1;
                l2*l1,m2*m1,n2*n1,m2*n1+m1*n2,    n2*l1+n1*l2,    l2*m1+l1*m2];

%T_stress=inv(T_strain)';
SM_rotated=T_strain*SM_origin*T_strain';

end
