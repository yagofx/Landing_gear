%--------------------------------------
%           CDIM Project
%--------------------------------------
clear all clc
format long



%Geomtrical data properties
r = 0.07367;    %m
r2 = 0.45;      %m
L1 = 0.63362;   %m
L2 = 0.40559;   %m
d1 = 0.388;     %m
d2 = 0.7344;    %m
d3 = 0.64;      %m
d4 = 0.57907;   %m
d5 = 0.54;      %m
d6 = 0.09;      %m

%Displacement function of s in time
smax = 0.36486; % m smax
beta = 50; %s time for smax
t_steps = 0.01;

t=0:t_steps:beta;
s = smax*(1/beta*t-1/(2*pi)*sin(2*pi/beta*t));
s_dot = smax/beta*(1-cos((2*pi/beta)*t));
s_ddot = 2*pi*smax/beta^2*sin((2*pi/beta)*t);

%Solving angular positions
q_vector = zeros(length(t),4); %generalized coordinates, phi1, phi2, phi3, phi4

for i=1:length(t)
    %Solution of the geometrical link equations
    pos = s(i);
    sol = @(q)landing_gear(q,pos,r,r2,L1,L2,d1,d2,d3,d4,d5,d6);

    if i == 1
        q0 =[1,1,1,1];
    else
        q0 = [q_vector(i-1,1),q_vector(i-1,2),q_vector(i-1,3),q_vector(i-1,4)];
    end
    q_vector(i,:) = fsolve(sol,q0);
end

%Separating angular positions
phi1 = q_vector(:,1);
phi2 = q_vector(:,2);
phi3 = q_vector(:,3);
phi4 = q_vector(:,4);

%Solving angular velocites
q_dot_d = zeros(length(t),4);

for i=1:length(t)
    phi_q_d = [sin(phi1(i))*(r + r2 + s(i)) d4*cos(phi2(i)) 0 0; -cos(phi1(i))*(r + r2 + s(i)) d4*sin(phi2(i)) 0 0; 0 L1*cos(phi2(i)) d5*cos(phi3(i)) -L2*sin(phi4(i)); 0 L1*sin(phi2(i)) -d5*sin(phi3(i)) L2*cos(phi4(i))];
    phi_t = [-cos(phi1(i))*s_dot(i); -sin(phi1(i))*s_dot(i); 0; 0];
    phi_q_i = [-cos(phi1(i)); -sin(phi1(i)); 0; 0];
    q_dot_i = s_dot(i);

    q_dot_d(i,:) = -inv(phi_q_d)*(phi_t+phi_q_i*q_dot_i);
end

%Separarating angular velocities
phi1_dot = q_dot_d(:,1);
phi2_dot = q_dot_d(:,2);
phi3_dot = q_dot_d(:,3);
phi4_dot = q_dot_d(:,4);


%Solving angular accelerations
q_ddot_d = zeros(length(t),4);

for i=1:length(t)
    phi_q_d = [sin(phi1(i))*(r + r2 + s(i)) d4*cos(phi2(i)) 0 0; -cos(phi1(i))*(r + r2 + s(i)) d4*sin(phi2(i)) 0 0; 0 L1*cos(phi2(i)) d5*cos(phi3(i)) -L2*sin(phi4(i)); 0 L1*sin(phi2(i)) -d5*sin(phi3(i)) L2*cos(phi4(i))];
    phi_q_i = [-cos(phi1(i)); -sin(phi1(i)); 0; 0];
    q_ddot_i = s_ddot(i);
    phi_q_dot = [0 sin(phi1(i))*s_dot(i) 0 0 0; 0 -cos(phi1(i))*s_dot(i) 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
    q_dot = [s_dot(i); phi1_dot(i); phi2_dot(i); phi3_dot(i); phi4_dot(i)];
    phi_t2 = [-cos(phi1(i))*s_ddot(i); -sin(phi1(i))*s_ddot(i); 0; 0];

    q_ddot_d(i,:) = -inv(phi_q_d)*(phi_q_i*q_ddot_i+phi_q_dot*q_dot+phi_t2);
end

%Separating angular accelerations
phi1_ddot = q_ddot_d(:,1);
phi2_ddot = q_ddot_d(:,2);
phi3_ddot = q_ddot_d(:,3);
phi4_ddot = q_ddot_d(:,4);

%Position plot
ax1 = nexttile;
plot(ax1,t,s)
title(ax1,'Displacement s')
ylabel(ax1,'$ s (m)$','Interpreter','latex')
xlabel(ax1,'$time (s)$','Interpreter','latex')

% Velocity plot
ax2 = nexttile;
plot(ax2,t,s_dot)
title(ax2,'Velocity s')
ylabel(ax2,'$ \dot{s} (m/s)$','Interpreter','latex')
xlabel(ax2,'$time (s)$','Interpreter','latex')
%Acceleration plot
ax3 = nexttile;
plot(ax3,t,s_ddot)
title(ax3,'Acceleration s')
ylabel(ax3,'$ \ddot{s} (m/s)$','Interpreter','latex')
xlabel(ax3,'$time (s)$','Interpreter','latex')

%Angular positions
ax4 = nexttile;
hold on
plot(ax4,t,phi1)
plot(ax4,t,phi2)
plot(ax4,t,phi3)
plot(ax4,t,phi4)
title(ax4,'Angular positions')
ylabel(ax4,'$ \phi (rad)$','Interpreter','latex')
xlabel(ax4,'$time (s)$','Interpreter','latex')
legend('\phi1','\phi2','\phi3','\phi4')
hold off

%Angular velocities
ax5 = nexttile;
hold on
plot(ax5,t,phi1_dot)
plot(ax5,t,phi2_dot)
plot(ax5,t,phi3_dot)
plot(ax5,t,phi4_dot)
title(ax5,'Angular velocities')
ylabel(ax5,'$\dot{\phi} (rad/s)$','Interpreter','latex')
xlabel(ax5,'$time (s)$','Interpreter','latex')
legend({'$\dot{\phi1}$','$\dot{\phi2}$','$\dot{\phi3}$','$\dot{\phi4}$'}, 'Interpreter', 'latex')
hold off

%Angular accelerations
ax6 = nexttile;
hold on
plot(ax6,t,phi1_ddot)
plot(ax6,t,phi2_ddot)
plot(ax6,t,phi3_ddot)
plot(ax6,t,phi4_ddot)
title(ax6,'Angular accelerations')
ylabel(ax6,'$\ddot{\phi} (rad/s^2)$','Interpreter','latex')
xlabel(ax6,'$time (s)$','Interpreter','latex')
legend({'$\ddot{\phi1}$','$\ddot{\phi2}$','$\ddot{\phi3}$','$\ddot{\phi4}$'}, 'Interpreter', 'latex')
hold off





