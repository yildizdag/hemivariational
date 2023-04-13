%---------------------------------------------------------------
% Single Spring
% Hemivariational Elasto-Damage Model
% Newton-Raphson Method
%---------------------------------------------------------------
k_el = 2;
kt = 1;
kd = 8;
F = 2.8;
f1 = @(u,d) k_el*(1-d)*u;
f2 = @(u,d) kt+kd*d-0.5*k_el*u^2;
%
u = 0;
d = 0;
del_u = 1;
del_d = 1;
TOL = 0.01;
count = 0;
while (abs(del_u)>TOL) && (abs(del_d)>TOL)
    tS = [k_el*(1-d), -k_el*u
             k_el*u, kd];
    del = tS\([F-f1(u,d); 0-f2(u,d)]);
    del_u = del(1);
    del_d = del(2);
    u = u + del_u;
    if del_d>0
        d = d + del_d;
    end
    count = count + 1;
end