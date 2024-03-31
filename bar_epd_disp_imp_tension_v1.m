%%%%%%%%%%%%%%%%%%%%%%%%%
% 1-D FEM
% Gaussian Quadrature
% Body force: f(x) = k^2 cos(2pikx/L)+3x
% ONLY Linear Elements
% Left end is FIXED (zero Dirichlet)
%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
% Parameters:
L = 1;
A = 1;
E = 1;
k = E*A;
kt = 1;
kd = 8;
st = 1.76;
sc = st;
ubar = 4;
Kp = 1E6*k;
ustep = 200;
% Quadrature Points:
ngp = 2;
[wgp, xgp] = gauss(ngp);
% Shape Functions and their derivs at quad. points:
N = [(1-xgp)/2 (1+xgp)/2];
dN = [-1/2*ones(length(xgp),1) 1/2*ones(length(xgp),1)];
% Number of elements:
ne = 1;
% Mesh (uniform elements):
nn = ne+1;
x = linspace(0,L,nn);
h = L/ne;
ndof = 3;
econn = zeros(ne,2);
sconn = zeros(ne,2*ndof);
dof = ndof*nn;
for e = 1:ne
    econn(e,:) = [e, e+1];
    sconn(e,:) = [6*e-5, 6*e-4, 6*e-3, 6*e-2, 6*e-1, 6*e];
end
%Storage
U = zeros(nn,ustep+1);
D = zeros(nn,ustep+1);
LT = zeros(nn,ustep+1);
LC = zeros(nn,ustep+1);
TOL = 1E-12;
count = 0;
for i=1:ustep
    % Displacement Step
    ui = (i/ustep)*ubar;
    U0 = U(:,i); U1 = U(:,i);
    D0 = D(:,i); D1 = D(:,i);
    LT0 = LT(:,i); LT1 = LT(:,i);
    LC0 = LC(:,i); LC1 = LC(:,i);
    enorm = 1;
    while enorm > TOL
        
        KT = zeros(dof,dof);
        R = zeros(dof,1);
        
        for e = 1:ne
            % Element Connectivity:
            lm = econn(e,:);
            slm = sconn(e,:);
            % Check State:
            if ( (1/h)*(U(lm(2),i)-U(lm(1),i))-0.5*(LT(lm(1),i)+LT(lm(2),i))+0.5*(LC(lm(1),i)+LC(lm(2),i)) ) >= 0
                
            elseif ( (1/h)*(U(lm(2),i)-U(lm(1),i))-0.5*(LT(lm(1),i)+LT(lm(2),i))+0.5*(LC(lm(1),i)+LC(lm(2),i)) ) < 0
                
            end
            
            if (((k/(2*kd))*(((1/h)*(U0(lm(2))-U0(lm(1)))-U0(lm(4))+U0(lm(5))))^2-kt/kd) < 0)
                ndof(lm(3)) = 0;
            end
            
            if ( ( ( (1/h)*(U0(lm(2))-U0(lm(1))) + U0(lm(5)) ) - st/(k*(1-U0(lm(3)))))  < 0) 
                ndof(lm(4)) = 0;
            end
            
            if ( ( -1*( (1/h)*(U0(lm(2))-U0(lm(1))) - U0(lm(4)) ) - sc/(k*(1-U0(lm(3)))))  < 0) 
                ndof(lm(5)) = 0;
            end
                
            for igp = 1:ngp
                Ng = N(igp,:);
                dNg = dN(igp,:);
                J = [x(lm(1)) x(lm(2))]*dNg';
                Nx = dNg./J;
                W = wgp(igp);
                Xg = [x(lm(1)) x(lm(2))]*Ng';
                kt_uu = (W*J).*((k*(1-U0(lm(3)))).*(Nx'*Nx));
                kt_ud = (W*J).*((-k*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))).*Nx');
                kt_ult = (W*J).*((-k*(1-U0(lm(3)))).*Nx');
                kt_ulc = (W*J).*((k*(1-U0(lm(3)))).*Nx');
                kt_dlt = (W*J)*(k*(Nx*[U0(lm(1));U0(lm(2))] - U0(lm(4)) + U0(lm(5))));
                kt_dlc = (W*J)*(-k*(Nx*[U0(lm(1));U0(lm(2))] - U0(lm(4)) + U0(lm(5))));
                kt_ltlc = (W*J)*(-k*(1-U0(lm(3))));
                
                KT([lm(1),lm(2)],[lm(1), lm(2)]) = KT([lm(1),lm(2)],[lm(1), lm(2)]) + kt_uu;
                KT([lm(1),lm(2)],lm(3)) = KT([lm(1),lm(2)],lm(3)) + kt_ud;
                KT(lm(3),[lm(1),lm(2)]) = KT(lm(3),[lm(1),lm(2)]) + kt_ud';
                KT([lm(1),lm(2)],lm(4)) = KT([lm(1),lm(2)],lm(4)) + kt_ult;
                KT(lm(4),[lm(1),lm(2)]) = KT(lm(4),[lm(1),lm(2)]) + kt_ult';
                KT([lm(1),lm(2)],lm(5)) = KT([lm(1),lm(2)],lm(5)) + kt_ulc;
                KT(lm(5),[lm(1),lm(2)]) = KT(lm(5),[lm(1),lm(2)]) + kt_ulc';
                KT(lm(3),lm(3)) = KT(lm(3),lm(3)) + (W*J)*kd;
                KT(lm(3),lm(4)) = KT(lm(3),lm(4)) + kt_dlt;
                KT(lm(4),lm(3)) = KT(lm(4),lm(3)) + kt_dlt;
                KT(lm(3),lm(5)) = KT(lm(3),lm(5)) + kt_dlc;
                KT(lm(5),lm(3)) = KT(lm(5),lm(3)) + kt_dlc;
                KT(lm(4),lm(4)) = KT(lm(4),lm(4)) + (W*J)*(k*(1-U0(lm(3))));
                KT(lm(4),lm(5)) = KT(lm(4),lm(5)) + kt_ltlc;
                KT(lm(5),lm(4)) = KT(lm(5),lm(4)) + kt_ltlc;
                KT(lm(5),lm(5)) = KT(lm(5),lm(5)) + (W*J)*(k*(1-U0(lm(3))));
                
                R([lm(1),lm(2)],1) = R([lm(1),lm(2)],1) + (W*J).*(((k*(1-U0(lm(3))))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))).*Nx');
                R(lm(3)) = R(lm(3)) + (W*J)*(-0.5*k*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))^2+kd*U0(lm(3))+kt);
                R(lm(4)) = R(lm(4)) + (W*J)*(-k*(1-U0(lm(3)))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))+st);
                R(lm(5)) = R(lm(5)) + (W*J)*(k*(1-U0(lm(3)))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))+sc);
                
            end
        end
        %Apply Dirichlet BCs:
        KT(1,1) = KT(1,1) + Kp;
        R(1) = R(1) + Kp*(U0(1));
        KT(nn,nn) = KT(nn,nn) + Kp;
        R(nn) = R(nn) + Kp*(U0(nn)-ui);
        %Delete Rows and Columns
        bc = find(ndof==0);
        act = setdiff(1:dof,bc);
        KT(bc,:) = [];
        KT(:,bc) = [];
        R(bc) = [];
        delU = KT\-R;
        for j = 1:length(act)
            if delU(j) > 0
                U1(act(j)) = U0(act(j)) + delU(j);
            end
        end
        e0 = bar_energy_epd(k,kd,kt,st,sc,Kp,ui,U0,x,conn);
        e1 = bar_energy_epd(k,kd,kt,st,sc,Kp,ui,U1,x,conn);
        enorm = abs(e1-e0)/e1;
        U0 = U1;
    end
    U(:,i+1) = U1;
end

% % % Loop over each element to assemble K and F:
% % for e = 1:ne
% %     lm = conn(e,:);
% %     for igp = 1:ngp
% %         Ng = N(igp,:);
% %         dNg = dN(igp,:);
% %         J = [x(lm(1)) x(lm(2))]*dNg';
% %         Nx = dNg./J;
% %         W = wgp(igp);
% %         Xg = [x(lm(1)) x(lm(2))]*Ng';
% %         ke = W*(E*A)*(Nx'*Nx)*J;
% %         fe = W*Ng'*f(Xg,k)*J;
% %         K(lm,lm) = K(lm,lm) + ke;
% %         F(lm) = F(lm) + fe;
% %     end
% % end
% % a = zeros(nn,1);
% % a(2:nn) = K(2:nn,2:nn)\F(2:nn);
% % %Plot:
% % figure;
% % plot(x,a)