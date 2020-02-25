% ADAMS_BASHFOURTH
%Adams-Bashforth mit k=1, p=1
a_bashforth_k1 = [-1, 1];
b_bashforth_k1 = [1, 0];

%Adams-Bashforth mit k=2, p=2
a_bashforth_k2 = [0, -1, 1];
b_bashforth_k2 = [-1/2, 3/2, 0];

%Adams-Bashforth mit k=3, p=3
a_bashforth_k3 = [0, 0, -1, 1];
b_bashforth_k3 = [5/12, -4/3, 23/12, 0];

%ADAMS_MOULTON
% Adams-Moulton mit k=1, p=2
a_moulton_k1 = [-1,1];
b_moulton_k1 = [1/2,1/2];

% Adams-Moulton mit k=2, p=3
a_moulton_k2 = [0, -1, 1];
b_moulton_k2 = [-1/12, 2/3, 5/12];

% Adams-Moulton mit k=3, p=4
a_moulton_k3 = [0, 0, -1, 1];
b_moulton_k3 = [1/24, -5/24, 19/24, 9/24];

%BDF
%Bdf mit k=1
a_bdf_k1 = [-1, 1]; 
b_bdf_k1 = [0, 1];

%Bdf mit k=2
a_bdf_k2 = [1/2, -2, 3/2]; 
b_bdf_k2 = [0, 0, 1];

%Bdf mit k=3
a_bdf_k3 = [-1/3, 3/2, -3, 11/6]; 
b_bdf_k3 = [0, 0, 0, 1];


