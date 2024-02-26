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
st = (1-0.1)*sqrt(2*k1*(0.1*kd+kt));
sc = st;
ubar = 4;
Kp = 1E6*k1;
N = 4000;
% Initial State (disp-damage-plastic/tension-plastic/comp)
u0 = [0, 0, 0, 0];
u1 = [0, 0, 0, 0];
% Storage
U = zeros(N+1,4);
F = zeros(N+1,1);
enorm = 1;
TOL = 1E-3;
count = 0;
for i=1:N
    % Displacement Step
    ui = (i/N)*ubar;
    while enorm > TOL
        if ((u0(1)-u0(3)+u0(4))>=0)
            KT = [k1*(1-u0(2))+Kp, -k1*(u0(1)-u0(3)+u0(4)), -k1*(1-u0(2));
                    -k1*(u0(1)-u0(3)+u0(4)), kd, k1*(u0(1)-u0(3)+u0(4));
                    -k1*(1-u0(2)), k1*(u0(1)-u0(3)+u0(4)), k1*(1-u0(2))];
            R = [k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+Kp*(u0(1)-ui);
                    -0.5*k1*(u0(1)-u0(3)+u0(4))^2+kd*u0(2)+kt;
                    -k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+st];
            dU = KT\-R;
            u1(1) = u0(1) + dU(1);
            if dU(2) > 0 && dU(3) <= 0
                u1(2) = u0(2) + dU(2);
                u1(3) = u0(3);
            end
            if dU(3) > 0 && dU(2) <= 0
                u1(3) = u0(3) + dU(3);
                u1(2) = u0(2);
            end
            if u1(2) > 1
                u1(2) = 1;
            end
            enorm = abs((e_epd(k1,kd,kt,st,sc,Kp,ui,u1)-e_epd(k1,kd,kt,st,sc,Kp,ui,u0))/e_epd(k1,kd,kt,st,sc,Kp,ui,u0));
            u0 = u1;
        elseif ((u0(1)-u0(3)+u0(4))<0)
            KT = [k2*(1-u0(2))+Kp, -k2*(u0(1)-u0(3)+u0(4)), k2*(1-u0(2))
                    -k2*(u0(1)-u0(3)+u0(4)), kd, -k2*(u0(1)-u0(3)+u0(4))
                      k2*(1-u0(2)), -k2*(u0(1)-u0(3)+u0(4)), k2*(1-u0(2))];
            R = [k2*(1-u0(2))*(u0(1)-u0(3)+u0(4))+Kp*(u0(1)-ui)
                    -0.5*k2*(u0(1)-u0(3)+u0(4))^2+kd*u0(2)+kt
                     k2*(1-u0(2))*(u0(1)-u0(3)+u0(4))+sc];
            dU = KT\-R;
            u1(1) = u0(1) + dU(1);
            if dU(2) > 0
                u1(2) = u0(2) + dU(2);
            else
                u1(2) = u0(2);
            end
            if dU(3) > 0
                u1(4) = u0(4) + dU(3);
            else
                u1(4) = u0(4);
            end
            u1(3) = u0(3);
            if u1(2) > 1
                u1(2) = 1;
            end
            enorm = abs((e_epd(k1,kd,kt,st,sc,Kp,ui,u1)-e_epd(k1,kd,kt,st,sc,Kp,ui,u0))/e_epd(k1,kd,kt,st,sc,Kp,ui,u0));
            u0 = u1;
        end
        count = count + 1;
    end
    U(i+1,:) = [u0(1), u0(2), u0(3), u0(4)];
    F(i+1) = k1*(1-u0(2))*(u0(1)-u0(3)+u0(4));
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