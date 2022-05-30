% GNU General Public License v3.0 | Copyright (c) 2022 Alexandros Pitilakis

clear all; close all; clc;

addpath( [pwd,'\Functions'] );

% Define an optical system. This can be either in ABCD form (+thickness) 
% but maximum functionality is offered when defining it in terms of optical
% surfaces (OS); each OS is a row in the allOS matrix where:
%   column #1: Absolute S-position of the OS
%   column #2: Optical power (1/Length) of OS. P=0 for plane surfaces.
%   column #3: Refr. index of material to the RIGHT of the OS
%   column #4: Transverse diameter of optical surface (lens). Rays hitting
%              the OS above/below the diameter will be stopped.
%  
% For example, a thick biconvex lens (R1=100mm, R2=50mm n=1.45, d=5mm) of
% diameter equal to 50mm in air is defined as:
%   OS(1,:) = [ 0 , (1.45-1)/+100 , 1.45 , 50 ];
%   OS(2,:) = [ 5 , (1-1.45)/-50  , 1.00 , 50 ];

% Define a single thin lens of focal length f=50, diameter=50, in air:
OS(1,:) = [ 0 , +1/50 , 1 , 50 ];
n12 = [ 1 1 ];

% Plot this optical system (with it's cardinal points etc) assuming it's
% working distance (WD=S1) is ~ 110
PlotParaxialOS( OS , n12 , 150 ) 

% Define a ray at distance S=-100, height=10, angle=-0.3rad, wavelength=633nm
rays = [ -100 , 10 , -0.3 , 633 ];

% Trace (until sfinal=+80) & plot this ray on the previously defined OS
sua = DoParaxialRayTracing( OS, n12, rays, 80 );

% Show the starting point of the ray with an asterisk
plot( rays(1,1) , rays(1,2) , 'r*' )






