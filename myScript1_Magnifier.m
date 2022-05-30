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
%
% ABCD matrices will be converted to their equivalent OS representation, 
% as thick lenses of appropriate {R1,R2,n} and a fixed diameter.
%
% Material indices (n1,n2) in the entrace and exit of OS are also required.
%
% You must also define a rectilinear (off-axis) object, with its S-position
% ATTENTION: S is in *absolute* NOT relative-to-A1 convention, meaning 
%            that negative S denotes a Real object, assuming that the
%            entrace of the optical system is at S=0.

% Choose a test case:
TestCase = 1;
switch TestCase
    case 1, % Thin diverging lens in air
        % Define ABCD
        ABCD = [ 1 0 ; +1/25 1];    n12 = [1 1]; d = 0; 
        % Define object S-position and u-height
        objS = +10; % virtual object!
        obju =   5;
        sfinal = [];
    case 2, % Exercises #3 (magnifier)
        ABCD = [0.6 3.8; -0.2 0.4]; n12 = [1 2]; d = 10; % Exercise #3
        objS = -1;
        obju =  5;
        sfinal = [];
end

% Convert ABCD matrix to optical surfaces (thick or thin lens)
if ~exist( 'OS' , 'var' )
    OS = Convert_ABCD2OS( ABCD , d , 0 );
end

% Call routine for image-formation using paraxially ray-traced cardinal
% rays (three of them)
ImageFormation_with_CardinalRayTracing( OS, n12, objS, obju, sfinal )


