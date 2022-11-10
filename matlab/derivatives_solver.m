%-------------------------------------
%           Derivatives solver
%-------------------------------------
clear all
format long

syms r r2 s(t) d1 d2 d3 d4 d5 d6 L1 L2 phi1 phi2 phi3 phi4 beta smax  pi t


f1 = -(r+r2+s(t))*cos(phi1)+d4*sin(phi2)-d6;
f2 = -(r+r2+s(t))*sin(phi1)-d4*cos(phi2)+d1;
f3 = -d2+d5*sin(phi3)+L2*cos(phi4)+L1*sin(phi2);
f4 = -d3+d5*cos(phi3)+L2*sin(phi4)-L1*cos(phi2)+d1;

%Geometrical link equations
phi= [f1, f2, f3, f4];

%Generalized coordinates
q = [s, phi1, phi2, phi3, phi4];

% Jacobian matrix
phi_q = jacobian(phi,q)

%Vector of temporal derivatives
phi_t = diff(phi,t)

%Temporal derivative of the jacobian
phi_q_dot = diff(phi_q,t)

%Vector of 2n temporal derivatives
phi_t2 = diff(phi_t,t)



