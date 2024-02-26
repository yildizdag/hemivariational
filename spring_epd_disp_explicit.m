%-------------------------------------------
%-------------------------------------------
% Single Spring Displacement Control
% Hemivariational Elasto-Damage Model
% Newton-Raphson Solver with Penalty Formulation
clear; clc;
%-------------------------------------------
%-------------------------------------------
k1=1;
k2=1;
kt = 1;
kd = 8;
st = 2.8;
sc = st;
ubar = 4;
N = 100;
% Initial State (disp-damage-plastic/tension-plastic/comp)
u = [0, 0, 0, 0];
% Storage
U = zeros(N+1,4);
F = zeros(N+1,1);
count = 0;
for i=1:N
    % Displacement Step
    ui = (i/N)*ubar;
    U(i+1,1) = ui;
    if (U(i+1,1)-U(i,3)+U(i,4)) >= 0
        if (((k1/(2*kd))*(U(i+1,1)-U(i,3)+U(i,4))^2-kt/kd)>0)
            U(i+1,2) = ((k1/(2*kd))*(U(i+1,1)-U(i,3)+U(i,4))^2-kt/kd);
        else
            U(i+1,2) = U(i,2);
        end
        if (U(i+1,1)+U(i,4)-st/(k1*(1-U(i+1,2))))>U(i,3)
            U(i+1,3) = U(i+1,1) + U(i,4) - st/(k1*(1-U(i+1,2)));
        else
            U(i+1,3) = U(i,3);
        end
        if (-(U(i+1,1)-U(i,3)) - sc/(k1*(1-U(i+1,2))))>U(i,4)
            U(i+1,4) = -(U(i+1,1)-U(i,3)) - sc/(k2*(1-U(i+1,2)));
        else
            U(i+1,4) = U(i,4);
        end
        F(i+1) = k1*(1-U(i+1,2))*(U(i+1,1)-U(i+1,3)+U(i+1,4));
    elseif (U(i+1,1)-U(i,3)+U(i,4)) < 0
        if (((k2/(2*kd))*(U(i+1,1)-U(i,3)+U(i,4))^2-(kt/kd))>=0)
            U(i+1,2) = ((k2/(2*kd))*(U(i+1,1)-U(i,3)+U(i,4))^2-kt/kd);
        else
            U(i+1,2) = U(i,2);
        end
        if (-(U(i+1,1)-U(i,3)) - sc/(k2*(1-U(i+1,2))))>U(i,4)
            U(i+1,4) = -(U(i+1,1)-U(i,3)) - sc/(k2*(1-U(i+1,2)));
        else
            U(i+1,4) = U(i,4);
        end
        if (U(i+1,1)+U(i,4)-st/(k2*(1-U(i+1,2))))>U(i,3)
            U(i+1,3) = U(i+1,1) + U(i,4) - st/(k2*(1-U(i+1,2)));
        else
            U(i+1,3) = U(i,3);
        end
        F(i+1) = k2*(1-U(i+1,2))*(U(i+1,1)-U(i+1,3)+U(i+1,4));
    end
    count = count + 1;
end
figure;
plot(U(:,1),U(:,2));
figure;
plot(U(:,1),U(:,3));
figure;
plot(U(:,1),U(:,4));
figure;
plot(U(:,1),F);