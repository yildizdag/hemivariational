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
ubar = 4;
Kpu = 1E9*k1;
Kpd = 1E9*k1;
Kpl = 1E9*k1;
N = 2000;
% Initial State (disp-damage-plastic/tension-plastic/comp)
u0 = [0, 0, 0, 0, 1E-6, 1E-6, 0];
% Storage
U = zeros(N+1,7);
F = zeros(N+1,1);
enorm = 1;
TOL = 1E-9;
count = 0;
for i=1:N
    % Displacement Step
    ui = (i/N)*ubar;
    di = U(i,2);
    li = U(i,3);
    while enorm > TOL
        KT = [k1*(1-u0(2))+Kpu, -k1*(u0(1)-u0(3)+u0(4)), -k1*(1-u0(2)), 0, 0;
                -k1*(u0(1)-u0(3)+u0(4)), Kpd+kd, k1*(u0(1)-u0(3)+u0(4)), -2*Kpd*u0(5), 0;
                -k1*(1-u0(2)), k1*(u0(1)-u0(3)+u0(4)), k1*(1-u0(2))+Kpl, 0, -2*Kpl*u0(6);
                  0, -2*u0(5)*Kpd, 0, 4*Kpd*u0(5)^2-2*Kpd*((u0(2)-di)-u0(5)^2), 0;
                  0, 0, -2*u0(6)*Kpl, 0, 4*Kpl*u0(6)^2-2*Kpl*((u0(3)-li)-u0(6)^2)];
        R = [k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+Kpu*(u0(1)-ui);
              -0.5*k1*(u0(1)-u0(3)+u0(4))^2+kd*u0(2)+kt+Kpd*((u0(2)-di)-u0(5)^2);
              -k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+st+Kpl*((u0(3)-li)-u0(6)^2);
              Kpd*((u0(2)-di)-u0(5)^2)*(-2*u0(5));
              Kpl*((u0(3)-li)-u0(6)^2)*(-2*u0(6))];
        dU = KT\-R;
        u0(1) = u0(1) + dU(1);
        u0(2) = u0(2) + dU(2);
        u0(3) = u0(3) + dU(3);
        u0(5) = u0(5) + dU(4);
        u0(6) = u0(6) + dU(5);
        enorm = norm(dU)
    end
    count = count + 1;
    U(i+1,:) = [u0(1), u0(2), u0(3), u0(4), u0(5), u0(6), u0(7)];
    F(i+1) = k1*(1-u0(2))*(u0(1)-u0(3)+u0(4));
    if U(i+1,5) < 1E-20
        u0(5) = 1E-10;
    end
    if U(i+1,6) < 1E-20
        u0(6) = 1E-10;
    end
    enorm = 1;
end
figure;
plot(U(:,1),U(:,2));
figure;
plot(U(:,1),U(:,3));
figure;
plot(U(:,1),F);