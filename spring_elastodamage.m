%---------------------------------------------------------------
% Single Spring
% Hemivariational Elasto-Damage Model
% Newton-Raphson Iterative Solver
% with Arc-Length Method
%---------------------------------------------------------------
k_el = 2;
kt = 1;
kd = 8;
F = 3;
R = 0.05;
lm0 = 0;
lm_max = 1;
u0 = 0;
u_max = 3;
f1 = @(u,d,lm) k_el*(1-d)*u-lm*F;
f2 = @(u,d,lm) kt+kd*d-0.5*k_el*u^2;
f3 = @(u,u0,d,lm,lm0) ((lm-lm0)^2)/(lm_max^2)+((u-u0)^2)/(u_max^2);

%
u = 0;
d = 0;
lm = 0;
del_u = 1;
del_d = 1;
del_lm = 1;
TOL = 0.001;
count = 0;
for i=1:step_max
    while (abs(del_u)>TOL) && (abs(del_d)>TOL && abs(del_lm)>TOL))
        tS = [k_el*(1-d), -k_el*u -F
            -k_el*u, kd, 0
            f3(u,u0,lm,lm0)^(-0.5)*((u-u0)/(u_max^2)), 0, f3(u,u0,lm,lm0)^(-0.5)*((lm-lm0)/(lm_max^2))];
        del = tS\([-f1(u,d,lm); -f2(u,d,lm); -(sqrt(f3(u,u0,d,lm,lm0))-R)]);
        del_u = del(1);
        del_d = del(2);
        del_lm = del(3);
        u = u + del_u;
        if del_d>0
            d = d + del_d;
        end
        lm = lm + del_lm;
        count = count + 1;
    end
    u0 = u;
    lm0 = lm;
end