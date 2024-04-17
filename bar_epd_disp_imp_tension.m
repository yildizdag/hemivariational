clear all; clc;
% Parameters:
L = 1;
A = 1;
E = 1;
k = E*A;
kt = 1;
kd = 8;
st = 1.2;
sc = st;
ubar = 0.2;
Kp = 1E9*k;
ustep = 2;
% Quadrature Points:
ngp = 2;
[wgp, xgp] = gauss(ngp);
% Shape Functions and their derivs at quad. points:
N = [(1-xgp)/2 (1+xgp)/2];
dN = [-1/2*ones(length(xgp),1) 1/2*ones(length(xgp),1)];
% Number of elements:
ne = 2;
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
TOL = 1E-12;
count = 0;
for i=1:ustep
    % Displacement Step
    ui = (i/ustep)*ubar;
    U0 = U(:,i); U1 = U(:,i);
    enorm = 1;
    while enorm > TOL
        
        KT = zeros(dof,dof);
        R = zeros(dof,1);
        ndof = 1:dof;
        
        for e = 1:ne
            lm = conn(e,:);
            if ( (1/h)*(U(lm(2),i)-U(lm(1),i))-U(lm(4),i)+U(lm(5),i) ) >= 0
                if ((((-k/2)*(((1/h)*(U(lm(2),i)-U(lm(1),i))-U(lm(4),i)+U(lm(5),i)))^2+kd*U(lm(3),i)+kt) > 0) && ((-k*(1-U(lm(3),i))*(((1/h)*(U(lm(2),i)-U(lm(1),i))-U(lm(4),i)+U(lm(5),i))) + st) > 0)) 
                    disp('damage and tension')
                    ndof(lm(5)) = 0;
                    
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
                        kt_dlt = (W*J)*(k*(Nx*[U0(lm(1));U0(lm(2))] - U0(lm(4)) + U0(lm(5))));
                        KT([lm(1),lm(2)],[lm(1), lm(2)]) = KT([lm(1),lm(2)],[lm(1), lm(2)]) + kt_uu;
                        KT([lm(1),lm(2)],lm(3)) = KT([lm(1),lm(2)],lm(3)) + kt_ud;
                        KT(lm(3),[lm(1),lm(2)]) = KT(lm(3),[lm(1),lm(2)]) + kt_ud';
                        KT([lm(1),lm(2)],lm(4)) = KT([lm(1),lm(2)],lm(4)) + kt_ult;
                        KT(lm(4),[lm(1),lm(2)]) = KT(lm(4),[lm(1),lm(2)]) + kt_ult';
                        KT(lm(3),lm(3)) = KT(lm(3),lm(3)) + (W*J)*kd;
                        KT(lm(3),lm(4)) = KT(lm(3),lm(4)) + kt_dlt;
                        KT(lm(4),lm(3)) = KT(lm(4),lm(3)) + kt_dlt;
                        KT(lm(4),lm(4)) = KT(lm(4),lm(4)) + (W*J)*(k*(1-U0(lm(3))));
                        R([lm(1),lm(2)],1) = R([lm(1),lm(2)],1) + (W*J).*(((k*(1-U0(lm(3))))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))).*Nx');
                        R(lm(3)) = R(lm(3)) + (W*J)*(-0.5*k*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))^2+kd*U0(lm(3))+kt);
                        R(lm(4)) = R(lm(4)) + (W*J)*(-k*(1-U0(lm(3)))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))+st);
                    end
                    
                elseif ((((-k/2)*(((1/h)*(U(lm(2),i)-U(lm(1),i))-U(lm(4),i)+U(lm(5),i)))^2+kd*U(lm(3),i)+kt) > 0) && ((-k*(1-U(lm(3),i))*(((1/h)*(U(lm(2),i)-U(lm(1),i))-U(lm(4),i)+U(lm(5),i))) + st) <= 0))
                    ndof(lm(4)) = 0; ndof(lm(5)) = 0;
                    disp('damage only')
                    for igp = 1:ngp
                        Ng = N(igp,:);
                        dNg = dN(igp,:);
                        J = [x(lm(1)) x(lm(2))]*dNg';
                        Nx = dNg./J;
                        W = wgp(igp);
                        Xg = [x(lm(1)) x(lm(2))]*Ng';
                        kt_uu = (W*J).*((k*(1-U0(lm(3)))).*(Nx'*Nx));
                        kt_ud = (W*J).*((-k*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))).*Nx');
                        KT([lm(1),lm(2)],[lm(1), lm(2)]) = KT([lm(1),lm(2)],[lm(1), lm(2)]) + kt_uu;
                        KT([lm(1),lm(2)],lm(3)) = KT([lm(1),lm(2)],lm(3)) + kt_ud;
                        KT(lm(3),[lm(1),lm(2)]) = KT(lm(3),[lm(1),lm(2)]) + kt_ud';
                        KT(lm(3),lm(3)) = KT(lm(3),lm(3)) + (W*J)*kd;
                        R([lm(1),lm(2)],1) = R([lm(1),lm(2)],1) + (W*J).*(((k*(1-U0(lm(3))))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))).*Nx');
                        R(lm(3),1) = R(lm(3),1) + (W*J)*(-0.5*k*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))^2+kd*U0(lm(3))+kt);  
                    end
                    
                elseif ((((-k/2)*(((1/h)*(U(lm(2),i)-U(lm(1),i))-U(lm(4),i)+U(lm(5),i)))^2+kd*U(lm(3),i)+kt) <= 0) && ((-k*(1-U(lm(3),i))*(((1/h)*(U(lm(2),i)-U(lm(1),i))-U(lm(4),i)+U(lm(5),i))) + st) > 0))
                    ndof(lm(3)) = 0; ndof(lm(5)) = 0;
                    disp('tension only')
                    for igp = 1:ngp
                        Ng = N(igp,:);
                        dNg = dN(igp,:);
                        J = [x(lm(1)) x(lm(2))]*dNg';
                        Nx = dNg./J;
                        W = wgp(igp);
                        Xg = [x(lm(1)) x(lm(2))]*Ng';
                        kt_uu = (W*J).*((k*(1-U0(lm(3)))).*(Nx'*Nx));
                        kt_ult = (W*J).*((-k*(1-U0(lm(3)))).*Nx');
                        KT([lm(1),lm(2)],[lm(1), lm(2)]) = KT([lm(1),lm(2)],[lm(1), lm(2)]) + kt_uu;
                        KT([lm(1),lm(2)],lm(4)) = KT([lm(1),lm(2)],lm(4)) + kt_ult;
                        KT(lm(4),[lm(1),lm(2)]) = KT(lm(4),[lm(1),lm(2)]) + kt_ult';
                        KT(lm(4),lm(4)) = KT(lm(4),lm(4)) + (W*J)*(k*(1-U0(lm(3))));                        
                        R([lm(1),lm(2)],1) = R([lm(1),lm(2)],1) + (W*J).*(((k*(1-U0(lm(3))))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))).*Nx');
                        R(lm(4)) = R(lm(4)) + (W*J)*(-k*(1-U0(lm(3)))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))+st);
                    end
                    
                else
                    ndof(lm(3)) = 0; ndof(lm(4)) = 0; ndof(lm(5)) = 0;
                    disp('elastic only')
                    for igp = 1:ngp
                        Ng = N(igp,:);
                        dNg = dN(igp,:);
                        J = [x(lm(1)) x(lm(2))]*dNg';
                        Nx = dNg./J;
                        W = wgp(igp);
                        Xg = [x(lm(1)) x(lm(2))]*Ng';
                        kt_uu = (W*J).*((k*(1-U0(lm(3)))).*(Nx'*Nx));
                        KT([lm(1),lm(2)],[lm(1), lm(2)]) = KT([lm(1),lm(2)],[lm(1), lm(2)]) + kt_uu;                      
                        R([lm(1),lm(2)],1) = R([lm(1),lm(2)],1) + (W*J).*(((k*(1-U0(lm(3))))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))).*Nx');
                    end
                end
             
            elseif ( ( (1/h)*(U(lm(2),i)-U(lm(1),i))-U(lm(4),i)+U(lm(5),i) ) < 0)
                
                if ((((-k/2)*(((1/h)*(U0(lm(2))-U0(lm(1)))-U0(lm(4))+U0(lm(5))))^2+kd*U0(lm(3))+kt) > 0) && ((k*(1-U0(lm(3)))*(((1/h)*(U0(lm(2))-U0(lm(1)))-U0(lm(4))+U0(lm(5)))) + sc) > 0) )
                    ndof(lm(4)) = 0;
                    for igp = 1:ngp
                        Ng = N(igp,:);
                        dNg = dN(igp,:);
                        J = [x(lm(1)) x(lm(2))]*dNg';
                        Nx = dNg./J;
                        W = wgp(igp);
                        Xg = [x(lm(1)) x(lm(2))]*Ng';
                        kt_uu = (W*J).*((k*(1-U0(lm(3)))).*(Nx'*Nx));
                        kt_ud = (W*J).*((-k*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))).*Nx');
                        kt_ulc = (W*J).*((k*(1-U0(lm(3)))).*Nx');
                        kt_dlc = (W*J)*(-k*(Nx*[U0(lm(1));U0(lm(2))] - U0(lm(4)) + U0(lm(5))));
                        KT([lm(1),lm(2)],[lm(1), lm(2)]) = KT([lm(1),lm(2)],[lm(1), lm(2)]) + kt_uu;
                        KT([lm(1),lm(2)],lm(3)) = KT([lm(1),lm(2)],lm(3)) + kt_ud;
                        KT(lm(3),[lm(1),lm(2)]) = KT(lm(3),[lm(1),lm(2)]) + kt_ud';
                        KT([lm(1),lm(2)],lm(5)) = KT([lm(1),lm(2)],lm(5)) + kt_ulc;
                        KT(lm(5),[lm(1),lm(2)]) = KT(lm(5),[lm(1),lm(2)]) + kt_ulc';
                        KT(lm(3),lm(3)) = KT(lm(3),lm(3)) + (W*J)*kd;
                        KT(lm(3),lm(5)) = KT(lm(3),lm(5)) + kt_dlc;
                        KT(lm(5),lm(3)) = KT(lm(5),lm(3)) + kt_dlc;
                        KT(lm(5),lm(5)) = KT(lm(5),lm(5)) + (W*J)*(k*(1-U0(lm(3))));
                        R([lm(1),lm(2)],1) = R([lm(1),lm(2)],1) + (W*J).*(((k*(1-U0(lm(3))))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))).*Nx');
                        R(lm(3)) = R(lm(3)) + (W*J)*(-0.5*k*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))^2+kd*U0(lm(3))+kt);
                        R(lm(5)) = R(lm(5)) + (W*J)*(k*(1-U0(lm(3)))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))+sc);
                    end
                    
                elseif ((((-k/2)*(((1/h)*(U0(lm(2))-U0(lm(1)))-U0(lm(4))+U0(lm(5))))^2+kd*U0(lm(3))+kt) > 0) && ((k*(1-U0(lm(3)))*(((1/h)*(U0(lm(2))-U0(lm(1)))-U0(lm(4))+U0(lm(5)))) + sc) <= 0) )
                    ndof(lm(4)) = 0; ndof(lm(5)) = 0;
                    for igp = 1:ngp
                        Ng = N(igp,:);
                        dNg = dN(igp,:);
                        J = [x(lm(1)) x(lm(2))]*dNg';
                        Nx = dNg./J;
                        W = wgp(igp);
                        Xg = [x(lm(1)) x(lm(2))]*Ng';
                        kt_uu = (W*J).*((k*(1-U0(lm(3)))).*(Nx'*Nx));
                        kt_ud = (W*J).*((-k*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))).*Nx');               
                        KT([lm(1),lm(2)],[lm(1), lm(2)]) = KT([lm(1),lm(2)],[lm(1), lm(2)]) + kt_uu;
                        KT([lm(1),lm(2)],lm(3)) = KT([lm(1),lm(2)],lm(3)) + kt_ud;
                        KT(lm(3),[lm(1),lm(2)]) = KT(lm(3),[lm(1),lm(2)]) + kt_ud';
                        KT(lm(3),lm(3)) = KT(lm(3),lm(3)) + (W*J)*kd;      
                        R([lm(1),lm(2)],1) = R([lm(1),lm(2)],1) + (W*J).*(((k*(1-U0(lm(3))))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))).*Nx');
                        R(lm(3)) = R(lm(3)) + (W*J)*(-0.5*k*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))^2+kd*U0(lm(3))+kt);              
                    end
                    
                elseif ((((-k/2)*(((1/h)*(U0(lm(2))-U0(lm(1)))-U0(lm(4))+U0(lm(5))))^2+kd*U0(lm(3))+kt) <= 0) && ((k*(1-U0(lm(3)))*(((1/h)*(U0(lm(2))-U0(lm(1)))-U0(lm(4))+U0(lm(5)))) + sc) > 0))
                    ndof(lm(3)) = 0; ndof(lm(4)) = 0;
                    for igp = 1:ngp
                        Ng = N(igp,:);
                        dNg = dN(igp,:);
                        J = [x(lm(1)) x(lm(2))]*dNg';
                        Nx = dNg./J;
                        W = wgp(igp);
                        Xg = [x(lm(1)) x(lm(2))]*Ng';
                        kt_uu = (W*J).*((k*(1-U0(lm(3)))).*(Nx'*Nx));
                        kt_ulc = (W*J).*((k*(1-U0(lm(3)))).*Nx');
                        KT([lm(1),lm(2)],[lm(1), lm(2)]) = KT([lm(1),lm(2)],[lm(1), lm(2)]) + kt_uu;
                        KT([lm(1),lm(2)],lm(3)) = KT([lm(1),lm(2)],lm(3)) + kt_ud;
                        KT(lm(3),[lm(1),lm(2)]) = KT(lm(3),[lm(1),lm(2)]) + kt_ud';
                        KT([lm(1),lm(2)],lm(5)) = KT([lm(1),lm(2)],lm(5)) + kt_ulc;
                        KT(lm(5),[lm(1),lm(2)]) = KT(lm(5),[lm(1),lm(2)]) + kt_ulc';
                        KT(lm(5),lm(5)) = KT(lm(5),lm(5)) + (W*J)*(k*(1-U0(lm(3))));
                        R([lm(1),lm(2)],1) = R([lm(1),lm(2)],1) + (W*J).*(((k*(1-U0(lm(3))))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))).*Nx');
                        R(lm(5)) = R(lm(5)) + (W*J)*(k*(1-U0(lm(3)))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))+sc);
                    end
                    
                else
                    ndof(lm(3)) = 0; ndof(lm(4)) = 0; ndof(lm(5)) = 0;
                    for igp = 1:ngp
                        Ng = N(igp,:);
                        dNg = dN(igp,:);
                        J = [x(lm(1)) x(lm(2))]*dNg';
                        Nx = dNg./J;
                        W = wgp(igp);
                        Xg = [x(lm(1)) x(lm(2))]*Ng';
                        kt_uu = (W*J).*((k*(1-U0(lm(3)))).*(Nx'*Nx));
                        KT([lm(1),lm(2)],[lm(1), lm(2)]) = KT([lm(1),lm(2)],[lm(1), lm(2)]) + kt_uu;
                        R([lm(1),lm(2)],1) = R([lm(1),lm(2)],1) + (W*J).*(((k*(1-U0(lm(3))))*(Nx*[U0(lm(1)); U0(lm(2))]-U0(lm(4))+U0(lm(5)))).*Nx');
                    end
                end
            end
        end
        %Apply Dirichlet BCs:
        KT(1,1) = KT(1,1) + Kp;
        R(1) = R(1) + Kp*(U0(1));
        KT(ne+1,ne+1) = KT(ne+1,ne+1) + Kp;
        R(ne+1) = R(ne+1) + Kp*(U0(ne+1)-ui);
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