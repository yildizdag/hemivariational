%---------------------------------------------------------------
%---------------------------------------------------------------
% Single Spring
% Hemivariational Elasto-Plastic-Damage Model
% Newton-Raphson Iterative Solver
%---------------------------------------------------------------
%---------------------------------------------------------------
k_el = 1;
kt = 1;
kd = 8;
dmax = 0.02;
st = (1-dmax)*sqrt(2*k_el*(dmax*kd+kt));
sc = st;
Kbar = 100000*k_el;
%
umax = 4;
N = 100;
ubar = [linspace(umax/(N),umax,N), linspace(umax-(umax/(N)),-umax,2*N+1), linspace(-umax+(umax/(N)),0,N)];
U = [0; 0; 0; 0];
delU = ones(4,1);
dStore = zeros(N,4);
TOL = 1E-3;
norm_dif = 1;
contatore = 0;
for i = 1:length(ubar)
    ubar_i = ubar(i);
    disp(ubar_i)
    while norm_dif > TOL
        if ((U(1)-U(3)+U(4))>=0)
            KT = [k_el*(1-U(2))+Kbar, -k_el*(U(1)-U(3)+U(4)), -k_el*(1-U(2)), k_el*(1-U(2));
                    -k_el*(U(1)-U(3)+U(4)), kd, k_el*(U(1)-U(3)+U(4)), -k_el*(U(1)-U(3)+U(4));
                    -k_el*(1-U(2)), k_el*(U(1)-U(3)+U(4)), k_el*(1-U(2)), -k_el*(1-U(2));
                    k_el*(1-U(2)), -k_el*(U(1)-U(3)+U(4)), -k_el*(1-U(2)), k_el*(1-U(2))];
            R = [k_el*(1-U(2))*(U(1)-U(3)+U(4))+Kbar*(U(1)-ubar_i); -0.5*k_el*(U(1)-U(3)+U(4))^2+kd*U(2)+kt; -k_el*(1-U(2))*(U(1)-U(3)+U(4))+st; k_el*(1-U(2))*(U(1)-U(3)+U(4))+sc];
            KT(end,:) = 0; KT(:,end) = 0; KT(end,end) = 1;
            R(end) = 0;
        else
            KT = [k_el*(1-U(2))+Kbar, -k_el*(U(1)-U(3)+U(4)), -k_el*(1-U(2)), k_el*(1-U(2));
                    -k_el*(U(1)-U(3)+U(4)), kd, k_el*(U(1)-U(3)+U(4)), -k_el*(U(1)-U(3)+U(4));
                    -k_el*(1-U(2)), k_el*(U(1)-U(3)+U(4)), k_el*(1-U(2)), -k_el*(1-U(2));
                      k_el*(1-U(2)), -k_el*(U(1)-U(3)+U(4)), -k_el*(1-U(2)), k_el*(1-U(2))];
            R = [k_el*(1-U(2))*(U(1)-U(3)+U(4))+Kbar*(U(1)-ubar_i); -0.5*k_el*(U(1)-U(3)+U(4))^2+kd*U(2)+kt; -k_el*(1-U(2))*(U(1)-U(3)+U(4))+st; k_el*(1-U(2))*(U(1)-U(3)+U(4))+sc];
            KT(3,:) = 0; KT(:,3) = 0; KT(3,3) = 1;
            R(3) = 0;
        end
        delU = KT\(-R);
        U(1) = U(1) + delU(1);
        if delU(2)>0
            U(2) = U(2) + delU(2);
        else
            delU(2) = 0;
        end
        if delU(3)>0
            U(3) = U(3) + delU(3);
        else
            delU(3) = 0;
        end
        if delU(4)>0
            U(4) = U(4) + delU(4);
        else
            delU(4) = 0;
        end
        norm_dif = norm(delU);
        contatore = contatore + 1;
    end
    delU = ones(4,1);
    norm_dif = 1;
    dStore(i+1,:) = [U(1), U(2), U(3), U(4)];
end
figure;
plot(dStore(:,1),k_el.*(1-dStore(:,2)).*(dStore(:,1)-dStore(:,3)+dStore(:,4)));
figure;
plot(dStore(:,1),dStore(:,2))
figure;
plot(dStore(:,1),dStore(:,3))