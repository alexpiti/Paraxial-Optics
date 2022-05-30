% GNU General Public License v3.0 | Copyright (c) 2022 Alexandros Pitilakis

clear all; close all; clc;

addpath( [pwd,'\Functions'] );

% -------------------------------------------------------------------------
% Keplerian Telescope
% -------------------------------------------------------------------------
n12 = [ 1 1 ];
OS(1,:) = [ 0   1/500 1   60 ]; % objective (fo=500mm)
OS(2,:) = [ 525 1/25  NaN 20 ]; % eyepiece (fe=25mm, dist 500+25mm)
sfinal = OS(end,1)+100;

% Plot the OS and annotations (cardinal points etc)
PlotParaxialOS( OS , n12 , 100 ); 

% Object is at infinity, so we define the so-called "telecentric" rays
% hitting the objective lens. We'll have to bundles, one //OA and one with
% a small tilt.
Nr = 10;
rays = ones(Nr,4);
rays(:,1) = -100;

% red rays are parallel to OA
rays(1:5,2) = OS(1,4)*linspace(-1/2, +1/2, Nr/2 );
rays(1:5,3) = 0;
rays(1:5,4) = 633;
% blue rays are coming from below the OA
rays(6:10,3) = 0.5*pi/180;
rays(6:10,2) = OS(1,4)*linspace(-1/2, +1/2, Nr/2 ) + rays(6:10,3)'.*rays(6:10,1)';
rays(6:10,4) = 480;

% Do ray tracing.
DoParaxialRayTracing( OS, n12, rays, sfinal);
axis normal;

title( 'Telescope (Kepler)')

% -------------------------------------------------------------------------
% Compound Microscope
% -------------------------------------------------------------------------
n12 = [ 1 1 ];
OS(1,:) = [ 0   1/20 1   10 ]; % objective (fo=20mm)
OS(2,:) = [ 150 1/30 NaN 10 ]; % eyepiece (fe=30mm, dist 150mm)
objS = -23.87; % position of object must be! accurately controlled!
obju = 0.1;
sfinal = OS(end,1)+30;

% Analytically calculate image position
ABCD = DoParaxialAnalysis( OS , n12 );
[S2,S1,mu] = DoParaxialImageFormation( ABCD , n12 , OS(end,1), 0, objS, [] );
if S2 > OS(end,1)
   error( 'Correct objS: No virtual image!' )
end

% Place object & image
figure;hold on; 
plot( objS*[1 1] , [0 obju   ] , 'ro-', 'MarkerFaceColor' , 'r' ); 
plot( S2  *[1 1] , [0 obju*mu] , 'go-', 'MarkerFaceColor' , 'g' )

% Plot the OS and annotations (cardinal points etc)
PlotParaxialOS( OS , n12 , 100 ); 

% Define Nr (e.g. 5) rays going into the object
Nr = 5;
rays = ones(Nr,4);
rays(:,1) = objS;
rays(:,2) = obju;
rays(:,3) = pi/20*linspace( -1/2, +1/2, Nr );
rays(:,4) = 633;

% Do ray tracing.
sua = DoParaxialRayTracing( OS, n12, rays, sfinal);
axis normal;

% Trace-back rays, from OS exit to virtual image (magnified)
for kr = 1 : Nr
    plot( [S2, sua(end-1,1,kr)] , [obju*mu, sua(end-1,2,kr)] , 'r:' )
end

title( 'Compound Microscope')
