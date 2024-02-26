%-------------------------------------------
%-------------------------------------------
% Single Spring Force Control
% Hemivariational Elasto-Plastic- Damage Model
% Fixed-Point Method
clear; clc;
%-------------------------------------------
%-------------------------------------------
k=1;
kt = 1;
kd = 8;
st = 1.8;
sc = st;
Fbar = 3;
% Initial State (disp-damage-plastic/tension-plastic/comp)
u0 = [0, 0, 0, 0, 0.001];
u1 = [0, 0, 0, 0, 0.001];
ubar = 0;
lmbar = 0;
R=0.005;
% Storage
step_max = 200;
U = zeros(step_max+1,5);
U(1,:) = u0;
%F = zeros(N+1,1);
enorm = 1;
TOL = 1E-4;
count = 0;
for i=1:step_max
    % Displacement Step
    %Fi = (i/N)*Fbar;
    while enorm > TOL
        u1(1) = (u0(5)*Fbar/(k*(1-u0(2))))+u0(3)-u0(4);
        u1(2) = (1/kd)*(0.5*k*(u0(1)-u0(3)+u0(4))^2-kt);
        u1(3) = (-st/(k*(1-u0(2))))+u0(1)+u0(4);
        u1(4) = (-sc/(k*(1-u0(2))))-u0(1)+u0(3);
        %u1(5) = lmbar + sqrt(R^2-(u0(1)-ubar)^2);
        u1(5) = (1/(2*lmbar))*(u0(5)^2+lmbar^2-R^2+(u0(1)-ubar)^2);
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
        % Energy Norm -- totalEnergy_epd(k,kd,kt,st,sc,Kp,ui,u0)
        totEn0 = totalEnergy_epd(k,kd,kt,st,sc,u0(5)*Fbar,u0);
        totEn1 = totalEnergy_epd(k,kd,kt,st,sc,u1(5)*Fbar,u1);
        enorm = abs(totEn1-totEn0)/totEn0;
        count = count + 1;
        u0 = u1;
    end
    U(i+1,:) = [u0(1), u0(2), u0(3), u0(4), u0(5)];
    %F(i+1) = k*(1-u0(2))*(u0(1)-u0(3)+u0(4));
    ubar = u0(1);
    lmbar = u0(5);
    u0(1) = u0(1) + 0.00001;
    u0(5) = u0(5) + 0.00001;
    enorm = 1;
end
figure;
plot(U(:,1),U(:,2));
figure;
plot(U(:,1),U(:,3));
figure;
plot(U(:,1),U(:,4));
%figure;
%plot(U(:,1),F);