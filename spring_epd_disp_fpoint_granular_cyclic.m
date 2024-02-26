%-------------------------------------------
%-------------------------------------------
% Single Spring Displacement Control
% Hemivariational Elasto-Damage Model
% Newton-Raphson Solver with Penalty Formulation
clear; clc;
%-------------------------------------------
%-------------------------------------------
kt=14E13;
kc=14E14;
Bt = 3.5E-7;
Bc = 1.5E-7;
st = 4E6;
sc = 4E7;
umax = 2E-7;
umin = -2E-7;
Kp = 1E9*kc;
N = 10;
ubar = [linspace(umax/(N),umax,N), linspace(umax-(umax/(N)),umin,2*N+1), linspace(umin+(umin/(N)),0,N)];
% Initial State (disp-damage-plastic/tension-plastic/comp)
u0 = [0, 0, 0, 0];
u1 = [0, 0, 0, 0];
% Storage
U = zeros(N+1,4);
F = zeros(N+1,1);
enorm = 1;
TOL = 1E-18;
count = 0;
for i=1:length(ubar)
    % Displacement Step
    ui = ubar(i);
    while enorm > TOL
        if (u0(1)-u0(3)+u0(4)) >= 0
            u1(1) = (1/(kt*(1-u0(2))+Kp))*(Kp*ui+(kt*(1-u0(2)))*(u0(3)-u0(4)));
            u1(2) = 1-exp(-(u0(1)-u0(3)+u0(4))/Bt);
        else
            u1(1) = (1/(kc*(1-u0(2))+Kp))*(Kp*ui+(kc*(1-u0(2)))*(u0(3)-u0(4)));
            u1(2) = (2/pi)*atan(-(u0(1)-u0(3)+u0(4))/Bc);
        end
        u1(3) = (-st/(kt*(1-u0(2))))+u0(1)+u0(4);
        u1(4) = (-sc/(kc*(1-u0(2))))-u0(1)+u0(3);
        if u1(2)<=u0(2)
            u1(2) = u0(2);
        end
        if u1(3)<=u0(3)
            u1(3) = u0(3);
        end
        if u1(4)<=u0(4)
            u1(4) = u0(4);
        end
        if u1(2) >= 1
            u1(2) = 0.99;
        end
        enorm = norm(u1-u0)
        count = count + 1;
        u0 = u1;
    end
    U(i+1,:) = [u0(1), u0(2), u0(3), u0(4)];
    if (u0(1)-u0(3)+u0(4)) >= 0
        F(i+1) = kt*(1-u0(2))*(u0(1)-u0(3)+u0(4));
    else
        F(i+1) = kc*(1-u0(2))*(u0(1)-u0(3)+u0(4));
    end
    enorm = 1;
end
figure;
plot(U(:,1),U(:,2));
figure;
plot(U(:,1),U(:,3));
figure;
plot(U(:,1),U(:,4));
figure;
plot(U(:,1),F);