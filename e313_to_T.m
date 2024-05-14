function [T] = e313_to_T(e313)
%E313_TO_T Convert Euler angles to attitude matrix using 313 sequence
% 
% Required
% --------
% e313 : double
%  (3,1) array of Euler angles given as [psi; theta; phi] in [rad]
%
% Returns
% -------
% T : double
%  (3,3) attitude matrix
%
% Written by Keith LeGrand, March 2023

psi = e313(1);
theta = e313(2);
phi = e313(3);

R_1 = @(a) [1,0,0;
             0,cos(a),sin(a);
             0,-sin(a),cos(a)];

R_3 = @(a) [cos(a),sin(a),0;
         -sin(a),cos(a),0;
         0,0,1];

T = R_3(phi)*R_1(theta)*R_3(psi);
end