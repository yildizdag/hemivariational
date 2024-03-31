function etot = bar_energy_epd(k,kd,kt,st,sc,Kp,ui,U0,x,conn)
% Quadrature Points:
ngp = 2;
[wgp, xgp] = gauss(ngp);
% Shape Functions and their derivs at quad. points:
dN = [-1/2*ones(length(xgp),1) 1/2*ones(length(xgp),1)];
etot = 0;
nn = size(conn,1)+1;
for e = 1:size(conn,1)
    lm = conn(e,:);
    for igp = 1:ngp
        dNg = dN(igp,:);
        J = [x(lm(1)) x(lm(2))]*dNg';
        Nx = dNg./J;
        W = wgp(igp);
        etot = etot + (W*J)*(0.5*k*(1-U0(lm(3)))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))^2+0.5*kd*U0(lm(3))^2+kt*U0(lm(3))+st*U0(lm(4))+sc*U0(lm(5)));
    end
end
etot = etot + 0.5*Kp*U0(1)^2 + 0.5*Kp*(U0(nn)-ui)^2;