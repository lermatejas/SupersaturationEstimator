function pid_sal = PID_Digital4(u_k)
% PID Digital (Astromand) 
% Eqs 8.23, 8.24 
K=-0.2;
Ti=10;
Td=10;
N=0.0015; % Filter Const. Derivative Action
b=1;      % Ref Weigth
ts=1;

%Reference
ysp=u_k(1);
%Process Variable
    y=u_k(2);
y_k_1=u_k(3);
% Error 
    e=u_k(4);
e_k_1=u_k(5);
% Action 
I_k_1=u_k(6)
D_k_1=u_k(7)

% Coefficients of different approximations of PID controller in continuous time
b_i1=0;
b_i2=(K*ts)/Ti;
a_d=1-((N*ts)/Td);
b_d=K*N;

%Control Law
P=K*(b*ysp-y);
I=I_k_1+b_i1*e+b_i2*e_k_1;
D=a_d*D_k_1-b_d*(y-y_k_1);
u_k=P+I+D;

% Saturation
if u_k<0
    u_k=0;
end
if u_k>0.015;
    u_k=0.015;
end
pid_sal=[u_k;I;D];