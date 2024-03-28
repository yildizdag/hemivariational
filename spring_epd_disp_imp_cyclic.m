%----------------------------------------------------
%----------------------------------------------------
% Single Spring Hemivariational Elasto-Plastic-Damage Model
% Displacement Control with Penalty Formulation 
% NEWTON-RAPHSON SOLVER
%----------------------------------------------------
clear; clc;
%----------------------------------------------------
k1 = 1;
k2 = 1;
kt = 1;
kd = 8;
st = 1.76;
sc = st;
ubar = 4;
Kp = 1E6*k1;
N = 500;
% Initial State (disp-damage-plastic/tension-plastic/comp)
u0 = [0, 0, 0, 0];
u1 = [0, 0, 0, 0];
% Storage
U = zeros(N+1,4);
F = zeros(N+1,1);
enorm = 1;
TOL = 1E-12;
count = 0;
for i=1:4*N
    % Displacement Step
    if i<N+1
        ui = (i/N)*ubar;
    elseif i<3*N+1
        ui = ubar-((i-N)/N)*ubar;
    elseif i<4*N+1
        ui = -ubar+((i-3*N)/N)*ubar;
    end
    while enorm > TOL
        % CHECK FIRST TENSION THEN KKT CONDITIONS
        if ((U(i,1)-U(i,3)+U(i,4))>=0)
            % BOTH DAMAGE AND PLASTICITY
            if ( (((((k1/(2*kd))*(u1(1)-u1(3)+u1(4))^2-kt/kd)))>0) && ((((u1(1)+u1(4)-st/(k1*(1-u1(2))))))>0) )
                 % TANGENT STIFFNESS
                 KT = [k1*(1-u0(2))+Kp, -k1*(u0(1)-u0(3)+u0(4)), -k1*(1-u0(2));
                         -k1*(u0(1)-u0(3)+u0(4)), kd, k1*(u0(1)-u0(3)+u0(4));
                         -k1*(1-u0(2)), k1*(u0(1)-u0(3)+u0(4)), k1*(1-u0(2))];
                 % RESIDUAL VECTOR
                 R = [k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+Kp*(u0(1)-ui);
                         -0.5*k1*(u0(1)-u0(3)+u0(4))^2+kd*u0(2)+kt;
                         -k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+st];
                 dU = KT\-R;
                 u1(1) = u0(1) + dU(1);
                 if dU(2) > 0
                    u1(2) = u0(2) + dU(2);
                 end
                 if dU(3) > 0
                    u1(3) = u0(3) + dU(3);
                 end
            % ONLY DAMAGE
            elseif ( (((((k1/(2*kd))*(u1(1)-u1(3)+u1(4))^2-kt/kd)))>0) && ((((u1(1)+u1(4)-st/(k1*(1-u1(2))))))<=0) )
                % TANGENT STIFFNESS
                KT = [k1*(1-u0(2))+Kp, -k1*(u0(1)-u0(3)+u0(4));
                        -k1*(u0(1)-u0(3)+u0(4)), kd];
                % RESIDUAL VECTOR
                R = [k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+Kp*(u0(1)-ui);
                    -0.5*k1*(u0(1)-u0(3)+u0(4))^2+kd*u0(2)+kt];
                dU = KT\-R;
                u1(1) = u0(1) + dU(1);
                if dU(2) > 0
                    u1(2) = u0(2) + dU(2);
                end
            % ONLY PLASTICITY
            elseif ( (((((k1/(2*kd))*(u1(1)-u1(3)+u1(4))^2-kt/kd)))<=0) && ((((u1(1)+u1(4)-st/(k1*(1-u1(2))))))>0) )
                %TANGENT STIFFNESS
                KT = [k1*(1-u0(2))+Kp, -k1*(1-u0(2));
                        -k1*(1-u0(2)), k1*(1-u0(2))];
                %RESIDUAL VECTOR
                R = [k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+Kp*(u0(1)-ui);
                      -k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+st];
                dU = KT\-R;
                u1(1) = u0(1) + dU(1);
                if dU(2) > 0
                    u1(3) = u0(3) + dU(2);
                end
            % ONLY ELASTIC
            else
                KT = k1*(1-u0(2))+Kp;
                R = k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+Kp*(u0(1)-ui);
                dU = -R/KT;
                u1(1) = u0(1) + dU;
            end
            enorm = abs((energy_epd(k1,kd,kt,st,sc,Kp,ui,u1)-energy_epd(k1,kd,kt,st,sc,Kp,ui,u0))/energy_epd(k1,kd,kt,st,sc,Kp,ui,u0));
            u0 = u1;
        % CHECK FIRST COMPRESSION THEN KKT CONDITIONS
        elseif ((U(i,1)-U(i,3)+U(i,4))<0)
            % BOTH DAMAGE AND PLASTICITY
            if ( (((((k1/(2*kd))*(u1(1)-u1(3)+u1(4))^2-kt/kd)))>0) && ( ((-(u1(1)-u1(3))-sc/(k1*(1-u1(2)))))>0) )
                 % TANGENT STIFFNESS
                 KT = [k1*(1-u0(2))+Kp, -k1*(u0(1)-u0(3)+u0(4)), k1*(1-u0(2));
                         -k1*(u0(1)-u0(3)+u0(4)), kd, -k1*(u0(1)-u0(3)+u0(4));
                           k1*(1-u0(2)), -k1*(u0(1)-u0(3)+u0(4)), k1*(1-u0(2))];
                 % RESIDUAL VECTOR
                 R = [k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+Kp*(u0(1)-ui);
                        -0.5*k1*(u0(1)-u0(3)+u0(4))^2+kd*u0(2)+kt;
                          k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+sc];
                 dU = KT\-R;
                 u1(1) = u0(1) + dU(1);
                 if dU(2) > 0
                    u1(2) = u0(2) + dU(2);
                 end
                 if dU(3) > 0
                    u1(4) = u0(4) + dU(3);
                 end
            % ONLY DAMAGE
            elseif ( (((((k1/(2*kd))*(u1(1)-u1(3)+u1(4))^2-kt/kd)))>0) && ( ((-(u1(1)-u1(3))-sc/(k1*(1-u1(2)))))<=0) )
                % TANGENT STIFFNESS
                KT = [k1*(1-u0(2))+Kp, -k1*(u0(1)-u0(3)+u0(4));
                        -k1*(u0(1)-u0(3)+u0(4)), kd];
                % RESIDUAL VECTOR
                R = [k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+Kp*(u0(1)-ui);
                    -0.5*k1*(u0(1)-u0(3)+u0(4))^2+kd*u0(2)+kt];
                dU = KT\-R;
                u1(1) = u0(1) + dU(1);
                if dU(2) > 0
                    u1(2) = u0(2) + dU(2);
                end
            % ONLY PLASTICITY
            elseif ( (((((k1/(2*kd))*(u1(1)-u1(3)+u1(4))^2-kt/kd)))<=0) && ( ((-(u1(1)-u1(3))-sc/(k1*(1-u1(2)))))>0) )
                %TANGENT STIFFNESS
                KT = [k1*(1-u0(2))+Kp, k1*(1-u0(2));
                          k1*(1-u0(2)), k1*(1-u0(2))];
                %RESIDUAL VECTOR
                R = [k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+Kp*(u0(1)-ui);
                        k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+sc];
                dU = KT\-R;
                u1(1) = u0(1) + dU(1);
                if dU(2) > 0
                    u1(4) = u0(4) + dU(2);
                end
            % ONLY ELASTIC
            else
                KT = k1*(1-u0(2))+Kp;
                R = k1*(1-u0(2))*(u0(1)-u0(3)+u0(4))+Kp*(u0(1)-ui);
                dU = -R/KT;
                u1(1) = u0(1) + dU;
            end
            enorm = abs((energy_epd(k1,kd,kt,st,sc,Kp,ui,u1)-energy_epd(k1,kd,kt,st,sc,Kp,ui,u0))/energy_epd(k1,kd,kt,st,sc,Kp,ui,u0));
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