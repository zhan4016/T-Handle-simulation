function [ax] = config_ax_for_3d(ax)
%CONFIG_AX_FOR_3D Configure plot axis for plotting rigid body
%
% Required
% --------
% ax : Axes
%   axis handle object for plotting on
%
% Returns
% -------
% ax : Axes
%  configured axis
%
% Written by Keith LeGrand, March 2023

axis(ax, 'equal')
axis(ax, [-10, 10, -10, 10, -10, 10])
ax.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
ax.YRuler.FirstCrossoverValue  = 0; % Y crossover with X axis
ax.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
ax.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
ax.XRuler.SecondCrossoverValue = 0; % X crossover with Z axis
ax.YRuler.SecondCrossoverValue = 0; % Y crossover with Z axis
view(ax, [105, 15])
axis(ax, 'off');
end