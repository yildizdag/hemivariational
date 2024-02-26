function etot = e_epd(k,kd,kt,st,sc,Kp,ui,u0)
etot = 0.5*k*(1-u0(2))*(u0(1)-u0(3)+u0(4))^2+0.5*Kp*(u0(1)-ui)^2+0.5*kd*u0(2)^2+kt*u0(2)+st*u0(3)+sc*u0(4);
end

