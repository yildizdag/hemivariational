%-------------------------------------------
%-------------------------------------------
% Single Spring Displacement Control
% Hemivariational Elasto-Damage Model
% Newton-Raphson Solver with Penalty Formulation
%-------------------------------------------
%-------------------------------------------
k = 1;
kt = 1;
kd = 8;
ubar = 4.0;
Kp = 1E6*k;
N = 100;
% Initial State (disp-damage)
u0 = [0, 0];
u1 = [0, 0];
% Storage
U = zeros(N+1,2);
F = zeros(N+1,1);
enorm = 1;
TOL = 1E-6;
count = 0;
for i=1:N
    % Displacement Step
    ui = (i/N)*ubar;
    while enorm > TOL
        KT = [k*(1-u0(2))+Kp, -k*u0(1)
                -k*u0(1),               kd];
        R = [k*(1-u0(2))*u0(1)+Kp*(u0(1)-ui);
                -0.5*k*u0(1)^2+kd*u0(2)+kt];
        dU = KT\-R;
        u1(1) = u0(1) + dU(1);
        if dU(2) > 0
            u1(2) = u0(2) + dU(2);
        end
        if u1(2) > 1
            u1(2) = 1;
        end
        enorm = abs((energy_elastodamage(k,kd,kt,Kp,ui,u1(1),u1(2))-energy_elastodamage(k,kd,kt,Kp,ui,u0(1),u0(2)))/energy_elastodamage(k,kd,kt,Kp,ui,u0(1),u0(2)));
        u0 = u1;
        count = count + 1;
    end
    U(i+1,:) = [u1(1), u1(2)];
    F(i+1) = k*(1-u1(2))*u1(1);
    enorm = 1;
end
figure;
plot(U(:,1),U(:,2));
figure;
plot(U(:,1),F);