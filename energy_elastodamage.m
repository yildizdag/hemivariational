function etot = energy_elastodamage(k,kd,kt,Kp,ui,u,d)
etot = 0.5*k*(1-d)*u^2+0.5*Kp*(u-ui)^2+0.5*kd*d^2+kt*d;
end

