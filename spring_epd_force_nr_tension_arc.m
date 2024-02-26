%-------------------------------------------
%-------------------------------------------
% Single Spring Force Control
% Hemivariational Elasto-Plastic- Damage Model
% Newton-Raphson Method
clear; clc;
%-------------------------------------------
%-------------------------------------------
k=1;
kt = 1;
kd = 8;
st = 1.8;
sc = st;
F = 3;
% Initial State (disp-damage-plastic/tension-plastic/comp)
u0 = [0, 0, 0, 0, 0.001];
u1 = [0, 0, 0, 0, 0];
S=0.1;
ubar = 0;
lbar = 0;
% Storage
step_max = 80;
umax = 4;
U = zeros(step_max+1,5);
%F = zeros(N+1,1);
enorm = 1;
TOL = 1E-3;
count = 0;
for i=1:step_max
    while enorm > TOL
        KT = [k*(1-u0(2)), -k*(u0(1)-u0(3)+u0(4)), -k*(1-u0(2)), -F;
                 -k*(u0(1)-u0(3)+u0(4)), kd, k*(u0(1)-u0(3)+u0(4)), 0;
                 -k*(1-u0(2)), k*(u0(1)-u0(3)+u0(4)), k*(1-u0(2)), 0;
                 (u0(1)-ubar)*((u0(1)-ubar)^2+(u0(5)-lbar)^2)^(-0.5), 0, 0, (u0(5)-lbar)*((u0(1)-ubar)^2+(u0(5)-lbar)^2)^(-0.5)];
        R = [k*(1-u0(2))*(u0(1)-u0(3)+u0(4))-u0(5)*F;
                -0.5*k*(u0(1)-u0(3)+u0(4))^2+kd*u0(2)+kt;
                -k*(1-u0(2))*(u0(1)-u0(3)+u0(4))+st;
                ((u0(1)-ubar)^2+(u0(5)-lbar)^2)^(0.5)-S];
        if (((k/(2*kd))*(u0(1)-u0(3)+u0(4))^2-(kt/kd)) < 0 && (-st/(k*(1-u0(2)))+u0(1)+u0(4)) < 0)
            KT(2:3,:) = []; KT(:,2:3) = [];
            R(2:3) = [];
            delU = KT\(-1.*R);
            u1(1) = u0(1) + delU(1);
            u1(5) = u0(5) + delU(2);
        elseif (((k/(2*kd))*(u0(1)-u0(3)+u0(4))^2-(kt/kd)) > 0 && (-st/(k*(1-u0(2)))+u0(1)+u0(4)) < 0)
            KT(3,:) = []; KT(:,3) = [];
            R(3) = [];
            delU = KT\(-1.*R);
            u1(1) = u0(1) + delU(1);
            u1(2) = u0(2) + delU(2);
            u1(5) = u0(5) + delU(3);
%             if delU(2) > 0
%                 u1(2) = u0(2) + delU(2);
%             else
%                 delU(2) = 0;
%             end
            if u1(2) >= 1
                u1(2) = 0.999;
            end
        elseif (((k/(2*kd))*(u0(1)-u0(3)+u0(4))^2-(kt/kd)) < 0 && (-st/(k*(1-u0(2)))+u0(1)+u0(4)) > 0)
            KT(2,:) = []; KT(:,2) = [];
            R(2) = [];
            delU = KT\(-1.*R);
            u1(1) = u0(1) + delU(1);
            u1(3) = u0(3) + delU(2);
            u1(5) = u0(5) + delU(3);
%             if delU(2) > 0
%                 u1(3) = u0(3) + delU(2);
%             else
%                 delU(2) = 0;
%             end
            if u1(2) >= 1
                u1(2) = 0.999;
            end
        else
            delU = KT\(-1.*R);
            u1(1) = u0(1) + delU(1);
            u1(2) = u0(2) + delU(2);
            u1(3) = u0(3) + delU(3);
            u1(5) = u0(5) + delU(4);
% %             if delU(2) > 0
% %                 u1(2) = u0(2) + delU(2);
% %             else
% %                 delU(2) = 0;
% %             end
% %             if delU(3) > 0
% %                 u1(3) = u0(3) + delU(3);
% %             else
% %                 delU(3) = 0;
% %             end
            if u1(2) >= 1
                u1(2) = 0.999;
            end
        end
        etot0 = totalEnergy_epd(k,kd,kt,st,sc,F,u0);
        etot1 = totalEnergy_epd(k,kd,kt,st,sc,F,u1);
        enorm = abs((etot0-etot1)/etot1);
        %enorm = norm(delU);
        count = count + 1;
        u0 = u1;
    end
    U(i+1,:) = [u0(1), u0(2), u0(3), u0(4), u0(5)];
    %F(i+1) = k*(1-u0(2))*(u0(1)-u0(3)+u0(4));
    ubar = u0(1);
    lbar = u0(5);
    u0(1) = u0(1) + 0.0000001;
    u0(5) = u0(5) + 0.0000001;
    enorm = 1;
    if U(i+1,1) > umax
        break;
    end
end
figure;
plot(U(1:i+1,1),U(1:i+1,2));
figure;
plot(U(1:i+1,1),U(1:i+1,3));
figure;
plot(U(1:i+1,1),F.*U(1:i+1,5));
%figure;
%plot(U(:,1),F);