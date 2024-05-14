function plot_T_handle(t,psi,theta,phi,student_name)
%PLOT_T_HANDLE Play animation of T handle given a time series of Euler
%angles.
%
% Required
% --------
% t : double
%  (nT,1) array of time values, where nT is the total number of time steps
% psi : double
%  (nT,1) array of first Euler rotation angle in [rad]
% theta : double
%  (nT,1) array of second Euler rotation angle in [rad]
% phi : double
%  (nT,1) array of third Euler rotation angle in [rad]
% student_name : char
%   first and last name of student, which will be added to the title of the
%   plot
%
% Written by Keith LeGrand and modified by Daniel Qi, March 2023

obj_import = stlread('T_handle.STL');
obj_fv.vertices = obj_import.Points/10; % in cm
obj_fv.faces = obj_import.ConnectivityList;
obj_fv = reducepatch(obj_fv,100);

% Move to (approx) CoM
obj_fv.vertices(:,1) = obj_fv.vertices(:,1) - 1;
obj_fv.vertices(:,3) = obj_fv.vertices(:,3) - 1;
obj_fv.vertices(:,2) = obj_fv.vertices(:,2) - 5;

% configure the figure
fig_inert_frame = figure;
% create a patch object of the object and configure appearance
ObjModel = patch(obj_fv);
ObjModel.FaceColor = [0.85, 0.33, 0.10];
ObjModel.FaceAlpha = 0.4;
% save the original coordinate values of the vertices
Vertices = ObjModel.Vertices;
ax_inert_frame = gca;
hold(ax_inert_frame,'on')
config_ax_for_3d(ax_inert_frame);
title(ax_inert_frame, student_name);
for k=1:length(t)
    T = e313_to_T([psi(k); theta(k); phi(k)]);
    % rotate vertices
    ObjModel.Vertices= (T.'*Vertices.').';
    pause(0.05)
end
end