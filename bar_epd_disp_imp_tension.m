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
kt = 0;
kd = 8;
st = 0.4;
sc = st;
ubar = 2;
Kp = 1E6*k;
ustep = 200;
% Quadrature Points:
ngp = 2;
[wgp, xgp] = gauss(ngp);
% Shape Functions and their derivs at quad. points:
N = [(1-xgp)/2 (1+xgp)/2];
dN = [-1/2*ones(length(xgp),1) 1/2*ones(length(xgp),1)];
% Number of elements:
ne = 20;
% Mesh (uniform elements):
nn = ne+1;
x = linspace(0,L,nn);
h = L/ne;
edof = 5;
conn = zeros(ne,edof);
dof = nn+3*ne;
for e = 1:ne
    conn(e,1:2) = [e, e+1];
    conn(e,3:5) = [(3*e-2)+nn, (3*e-1)+nn, (3*e)+nn];
end
%Storage
U = zeros(dof,ustep+1);
F = zeros(ustep+1,1);
TOL = 1E-6;
count = 0;
for i=1:ustep
    % Displacement Step
    ui = (i/ustep)*ubar;
    U0 = U(:,i);
    enorm = 1;
    while enorm > TOL
        KT = zeros(dof,dof);
        R = zeros(dof,1);
        ndof = 1:dof;
        for e = 1:ne
            
            lm = conn(e,:);
            
            if (((0.5*k*(((1/h)*(U0(lm(2))-U0(lm(1)))-U0(lm(4))+U0(lm(5))))^2-kd*U0(lm(3))-kt)) < 0)
                ndof(lm(3)) = 0;
            end
            
            if ( ( (k*(1-U0(lm(3))))*( (1/h)*(U0(lm(2))-U0(lm(1))) - U0(lm(4)) + U0(lm(5)) ) - st)  < 0) 
                ndof(lm(4)) = 0;
            end
            
            if ( ( (-k*(1-U0(lm(3))))*( (1/h)*(U0(lm(2))-U0(lm(1))) - U0(lm(4)) + U0(lm(5)) ) - sc)  < 0) 
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
                U0(act(j)) = U0(act(j)) + delU(j);
            end
        end
        enorm = norm(delU)
    end
    U(:,i+1) = U0;
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