%% AAE340 HW12 Main
% Michael Zhang

clear; close all;
%% Initi
% Simulation Parameters
max_time = 10; %[sec] simulation time
int_time = 0.07; %[sec] time interval
t_lin = 0:int_time:max_time;

psi_0 = 0; %[rad]
theta_0 = pi/2; %[rad]
phi_0 = 0; %[rad]

%major axis                    minor axis
omega1_0 = 0.05; omega2_0 = 8; omega3_0 = 0.1; %[rad/sec] (unstable) <<<<
%omega1_0 = 0.1; omega2_0 = 0.05; omega3_0 = 8; %[rad/sec] (stable) <<<<

options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);

% T-Handle
C1_H = 10;   %[cm] c1 height
C1_D = 2;    %[cm] c1 diameter
C1_M = 50;   %[g] c1 mass

C2_H = 2;    %[cm] c2 height
C2_D = 1;    %[cm] c2 diameter
C2_M = 20;   %[g] c2 mass 

% center of mass
y_cg = (C1_M*(C2_H+C1_D/2)+C2_M*C2_H/2)/(C1_M+C2_M); %[cm]
x_cg = 0; %[cm]
z_cg = 0; %[cm]
mass_tot = C1_M+C2_M; %[g]

% moment of inertia
%C2
Iy2 = 0.5*C2_M*(C2_D/2)^2;
Ix2 = 0.25*C2_M*(C2_D/2)^2+1/12*C2_M*C2_H^2;
Iz2 = Ix2;

%C1
Iy1 = 0.25*C1_M*(C1_D/2)^2+1/12*C1_M*C1_H^2;
Ix1 = Iy1;
Iz1 = 0.5*C1_M*(C1_D/2)^2;

% Total MOI
r_c1_cg = C2_H+C1_D/2-y_cg;
r_c2_cg = y_cg-C2_H/2;
Ix = r_c1_cg^2*C1_M+Ix1 + r_c2_cg^2*C2_M+Ix2; % Parallel Axis Theorem
Iy = Iy1+Iy2;
Iz = r_c1_cg^2*C1_M+Iz1 + r_c2_cg^2*C2_M+Iz2;

I = [Ix 0 0; 0 Iy 0; 0 0 Iz]; % inertia matrix

%display
fprintf("Moment of inertia of cylinder 1 is given by: \nIx = %.3f [g-cm^2] \nIy = %.3f [g-cm^2] \nIz = %.3f [g-cm^2]\n",Ix1,Iy1,Iz1)
fprintf("Moment of inertia of cylinder 2 is given by: \nIx = %.3f [g-cm^2] \nIy = %.3f [g-cm^2] \nIz = %.3f [g-cm^2]\n",Ix2,Iy2,Iz2)
fprintf("The distance betweent the C/M of Cylinder 1 and the C/M of the T-Handle is %.3f [cm]\n",r_c1_cg)
fprintf("The distance betweent the C/M of Cylinder 2 and the C/M of the T-Handle is %.3f [cm]\n",r_c2_cg)
fprintf("The inertia tensor of the T-Handle is:\n")
disp(I)


%% Simulation
eq1 = @(t,x) [1/sin(x(2))*(x(4)*sin(x(3))+x(5)*cos(x(3))); ...
    x(4)*cos(x(3))-x(5)*sin(x(3)); ...
    x(6)-cos(x(2))/sin(x(2))*(x(4)*sin(x(3))+x(5)*cos(x(3)));...
    -(I(3,3)-I(2,2))*x(5)*x(6)/I(1,1);...
    -(I(1,1)-I(3,3))*x(6)*x(4)/I(2,2);...
    -(I(2,2)-I(1,1))*x(5)*x(4)/I(3,3)];
[t, ca1_nu] = ode45(eq1,t_lin,[phi_0 theta_0 psi_0 omega1_0 omega2_0 omega3_0], options);

fprintf("The inital angular velocities are:\n omega1 = %.f [rad/sec]\n omega2 = %.f [rad/sec]\n omega3 = %.f [rad/sec]\n",omega1_0,omega2_0,omega3_0);

phi = ca1_nu(:,1); theta = ca1_nu(:,2); psi = ca1_nu(:,3); %[rad]
omega1 = ca1_nu(:,4); omega2 = ca1_nu(:,5); omega3 = ca1_nu(:,6); %[rad/sec]

%% Vector calculation
omega_B = [omega1 omega2 omega3]; %intial angular velocity vector
H_B = [omega1*I(1,1) omega2*I(2,2) omega3*I(3,3)]; %intial angular velocity vector

%preallocate
omega_I = zeros([length(psi) 3]); %angular velocity B frame
H_I = zeros([length(psi) 3]); %angular momentum B frame

%transformation
for i = 1:length(psi)
    l = e313_to_T([phi(i) theta(i) psi(i)]);
    omega_I(i,:) = inv(l)*omega_B(i,:)';
    H_I(i,:) = inv(l)*H_B(i,:)';
end

figure(2)
plot(t,omega_B)
title({"(c).i ^I\omega^{B} in B frame"},{"-Michael Zhang"});grid minor;
ylabel('^I\omega^{B} [rad/sec]');xlabel('Time [sec]')
legend(["\omega_1", "\omega_2", "\omega_3"])

figure(3)
plot(t,omega_I)
title({"(c).ii ^I\omega^{B} in I frame"},{"-Michael Zhang"});grid minor;
ylabel('^I\omega^{B} [rad/sec]');xlabel('Time [sec]')
legend(["H_1", "H_2", "H_3"])

figure(4)
plot(t,H_B)
title({"(c).iii ^IH^{B} in B frame"},{"-Michael Zhang"});grid minor;
ylabel('^IH^{B} [g cm^2/sec]');xlabel('Time [sec]')
legend(["\omega_1", "\omega_2", "\omega_3"])

figure(5)
plot(t,H_I)
title({"(c).iv ^IH^{B} in I frame"},{"-Michael Zhang"});grid on;
ylabel('^IH^{B} [g cm^2/sec]');xlabel('Time [sec]')
legend(["H_1", "H_2", "H_3"])

plot_T_handle(t,ca1_nu(:,3),ca1_nu(:,2),ca1_nu(:,1),'T-Handle Simulation - Michael Zhang');



































