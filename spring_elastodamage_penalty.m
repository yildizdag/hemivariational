%---------------------------------------------------------------
%---------------------------------------------------------------
% Single Spring
% Hemivariational Elasto-Damage Model
% Newton-Raphson Iterative Solver
% Penalty Formulation for Irreversibility
clear; clc;
%---------------------------------------------------------------
%---------------------------------------------------------------
k_el = 1;
kt = 1;
kd = 8;
K1 = 1000*k_el;
K2 = 100*kd;
ubar = 2.0;
N = 80;
U = [0; 0; 1];
delU = ones(3,1);
dStore = zeros(N,3);
TOL = 1E-12;
norm_dif = 1;
contatore = 0;
for i=1:N
    ubar_i = (i/N)*ubar;
    d0 = U(2);
    disp(i);
    while (norm(delU(1:2))>TOL)
        KT = [k_el*(1-U(2))+K1, -k_el*U(1), 0; 
                -k_el*U(1), kd+K2, -2*K2*U(3);
                0, -2*K2*U(3), -2*K2*(U(2)-d0)+6*K2*U(3)^2];
        R = [k_el*(1-U(2))*U(1)+K1*(U(1)-ubar_i); -0.5*k_el*U(1)^2+kd*U(2)+kt+K2*((U(2)-d0)-U(3)^2); -2*K2*U(3)*((U(2)-d0)-U(3)^2)];
        delU = KT\(-R);
        U = U + delU;
        if U(2) >= 1.0
            U(2) = 1.0;
        end
    end
    delU = ones(3,1);
    dStore(i,:) = [U(1), U(2), U(3)];
    contatore = contatore + 1;
end
figure;
plot([0;dStore(:,1)],[0;dStore(:,2)])
figure;
plot([0;dStore(:,1)],[0;dStore(:,3)])
disp(contatore);