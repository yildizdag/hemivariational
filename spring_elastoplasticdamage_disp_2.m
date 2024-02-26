%-------------------------------------------
%-------------------------------------------
% Single Spring Displacement Control
% Hemivariational Elasto-Damage Model
% Newton-Raphson Solver with Penalty Formulation
clear; clc;
%-------------------------------------------
%-------------------------------------------
k1 = 1;
k2 = 1;
kt = 1;
kd = 8;
st = (1-0.2)*sqrt(2*k1*(0.2*kd+kt));
sc = st;
ubar = 0.1;
Kp = 1E6*k1;
N = 10;
% Initial State (disp-damage-plastic/tension-plastic/comp)
u0 = [0, 0, 0, 0];
u1 = [0, 0, 0, 0];
% Storage
U = zeros(N+1,4);
F = zeros(N+1,1);
enorm = 1;
TOL = 1E-8;
count = 0;
for i=1:N
    % Displacement Step
    ui = (i/N)*ubar;
    while enorm > TOL
%         if ((u0(1)-u0(3)+u0(4))>=0)
        KT = [k1*(1-u0(2))+Kp, -k1*(u0(1)-u0(3)+u0(4)), -k1*(1-u0(2)), k2*(1-u0(2));
                    -k1*(u0(1)-u0(3)+u0(4)), kd, k1*(u0(1)-u0(3)+u0(4)), -k2*(u0(1)-u0(3)+u0(4));
                    -k1*(1-u0(2)), k1*(u0(1)-u0(3)+u0(4)), k1*(1-u0(2)), -k2*(1-u0(2));
                    -k1*(1-u0(2)), -k1*(u0(1)-u0(3)+u0(4)), -k1*(1-u0(2)), k2*(1-u0(2))];
        R = [k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+Kp*(u0(1)-ui)
            -0.5*k1*(u0(1)-u0(3)+u0(4))^2+kd*u0(2)+kt
            -k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+st;
                k2*(1-u0(2))*(u0(1)-u0(3)+u0(4))+sc];
        dU = KT\-R;
        u1 = u0 + dU';
        enorm = norm(dU);
        u0 = u1;
        count = count + 1;
    end
    U(i+1,:) = [u1(1), u1(2), u1(3), u1(4)];
    enorm = 1;
end