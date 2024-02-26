%-------------------------------------------
%-------------------------------------------
% Single Spring Displacement Control
% Hemivariational Elasto-Damage Model
% Newton-Raphson Solver with Penalty Formulation
clear; clc;
%-------------------------------------------
%-------------------------------------------
k = 1;
kt = 1;
kd = 8;
ubar = 4.0;
Kpu = 1E9*k;
Kpd = 1E9*kt;
N = 200;
% Initial State (disp-damage)
u0 = [0, 0, 1E-6];
u1 = [0, 0, 0];
% Storage
U = zeros(N+1,3);
F = zeros(N+1,1);
enorm = 1;
TOL = 1E-8;
count = 0;
for i=1:N
    % Displacement Step
    ui = (i/N)*ubar;
    di = U(i,2);
    while enorm > TOL
        KT = [k*(1-u0(2))+Kpu, -k*u0(1), 0;
                -k*u0(1),               kd+Kpd, -2*Kpd*u0(3);
                0, -2*u0(3)*Kpd, 4*Kpd*u0(3)^2-2*Kpd*((u0(2)-di)-u0(3)^2)];
        R = [k*(1-u0(2))*u0(1)+Kpu*(u0(1)-ui);
              -0.5*k*u0(1)^2+kd*u0(2)+kt+Kpd*((u0(2)-di)-u0(3)^2);
                Kpd*((u0(2)-di)-u0(3)^2)*(-2*u0(3))];
        dU = KT\-R;
        u1 = u0 + dU';
        enorm = norm(dU);
        %enorm = abs((energy_elastodamage(k,kd,kt,Kp,ui,u1(1),u1(2))-energy_elastodamage(k,kd,kt,Kp,ui,u0(1),u0(2)))/energy_elastodamage(k,kd,kt,Kp,ui,u0(1),u0(2)));
        u0 = u1;
        count = count + 1;
    end
    U(i+1,:) = [u1(1), u1(2), u1(3)];
    F(i+1) = k*(1-u1(2))*u1(1);
    enorm = 1;
end
figure;
plot(U(:,1),U(:,2));
figure;
plot(U(:,1),F);