% GNU General Public License v3.0 | Copyright (c) 2022 Alexandros Pitilakis

clear all; close all; clc;

addpath( [pwd,'\Functions'] );

%% === Optical System ===

% Define a 4-lens 2-stop OS: [S-pos, OptPower, refr.indx, diameter]
OS(1,:) = [  0 , +1/20 , 1.00 , 30 ]; % thin lens: f=1/P=+20, diameter=30
OS(2,:) = [  4 , -1/15 , 1.00 , 20 ];
OS(3,:) = [ 10 , 0     , 1.00 , 15 ]; % stop #1 (AS)
OS(4,:) = [ 16 , +1/15 , 1.00 , 20 ];
OS(5,:) = [ 20 , +1/40 , 1.00 , 30 ];
OS(6,:) = [ 31 , 0     , 1.00 , 10 ]; % stop #2 (FS)
n12 = [ 1 1 ];
WD = 100; % working-distance "S1" (positive == to the LEFT of first lens)
sfinal = 50; % final S-position (absolute) for ray plots

% Paraxial analysis of this optical system returns: its total ABCD matrix,
% its cardinal points (CPs), its Pupils/AS and Windows/FS.
[ABCD,CPs,Pupils,Windows] = DoParaxialAnalysis( OS , n12 , WD );



%% === Simple object formation ===
figure;

% Plot OS etc
PPOSsettings.PlotCPs = 1; % plot cardinal points (F,H,N)
PPOSsettings.PlotApertureAnalysis = [ 0 0 ]; % Pupils,Windows
PPOSsettings.CPs_Legend = 1; % Plot cardinal points (F,H,N) legend
PPOSsettings.ApA_Legend = [0 0]; % Pupils+Windows legend
PlotParaxialOS( OS , n12 , WD , PPOSsettings ); 

% Define a bunch of rays from a point object to the first-lens aperture
obju = -10;
Nr = 11;
rays = NaN*ones(Nr,4);
rays(:,1) = -WD; % S-position (axial)
rays(:,2) =  obju; % u-position (height)
angtop = ( +OS(1,4)/2 - obju ) / WD;
angbot = ( -OS(1,4)/2 - obju ) / WD;
rays(:,3) = linspace( angbot , angtop, Nr ); % angle [rad]
rays(:,4) = round( linspace( 450,633,Nr) ); % wavelength (color) [nm]

% Do ray-tracing and plot rays
sua = DoParaxialRayTracing( OS, n12, rays, sfinal, 0 , 1);
%pause(30)
for kr = 1 : Nr
    plot( sua(:,1,kr), sua(:,2,kr), '-' , 'Color' , ...
        calculateVisibleSpectrumColor(rays(kr,4))  );
end

% Now, define rays coming from object@ infinity and hitting the first-lens
Nr = 11;
rays = NaN*ones(Nr,4);
rays(:,1) = -WD; % S-position (axial)
rays(:,3) = -10*pi/180; % angle [rad]
rays(:,2) = rays(:,1).*rays(:,3) + OS(1,4)*linspace(-1/2,+1/2,Nr)'; % u-position (height)
rays(:,4) = 700; % wavelength (color) [nm]

% Do ray-tracing and plot rays
sua = DoParaxialRayTracing( OS, n12, rays, sfinal, 0 , 1);
%pause
for kr = 1 : Nr
    plot( sua(:,1,kr), sua(:,2,kr), '-' , 'Color' , ...
        calculateVisibleSpectrumColor(rays(kr,4))  );
end

title( 'Image Formation' );



%% === Pupils (and aperture-stop, AS) analysis ===
figure;

% Plot OS etc
PPOSsettings.PlotCPs = 0; % plot cardinal points (F,H,N)
PPOSsettings.PlotApertureAnalysis = [ 1 0 ]; % Pupils,Windows
PPOSsettings.CPs_Legend = 0; % Plot cardinal points (F,H,N) legend
PPOSsettings.ApA_Legend = [1 0]; % Pupils+Windows legend
PlotParaxialOS( OS , n12 , WD , PPOSsettings ); 

% Define a bunch of rays from the axial-point of object plane, 
% with increasing ray-angle (tilt)
Nr = 20;
rays = NaN*ones(Nr,4);
rays(:,1) = -WD; % S-position (axial)
rays(:,2) =   0; % u-position (height)
rays(:,3) = 30 * linspace(-1/2,1/2,Nr)' * pi/180; % angle [rad]
wls = round( linspace( 633, 455, Nr/2) )'; % wavelength (color) [nm]
rays(:,4) = [wls; flipud(wls)];

% Do ray-tracing and plot rays
sua = DoParaxialRayTracing( OS, n12, rays, sfinal, 0 , 1);
for kr = 1 : Nr
    plot( sua(:,1,kr), sua(:,2,kr), '-' , 'Color' , ...
        calculateVisibleSpectrumColor(rays(kr,4))  );
end

% Define and plot limiting rays, from object-plane axial-point to P_en edges
raysAS = NaN*ones(2,4);
raysAS(:,1) = -WD; % S-position (axial)
raysAS(:,2) =   0; % u-position (height)
raysAS(:,3) = Pupils.enD/(+WD+Pupils.enS) * [-1/2 +1/2]; % angle [rad]
raysAS(:,4) = 700; % wavelength (color) [nm]
suaAS = DoParaxialRayTracing( OS, n12, raysAS, sfinal, 0 , 1);
for kr = 1 : size(raysAS,1)
    plot( suaAS(:,1,kr), suaAS(:,2,kr), '--' , 'Color' , ...
        calculateVisibleSpectrumColor(raysAS(kr,4))  );
end

title( 'Pupils / Aperture-Stop analysis' );



%% === Windows (and Field-Stop, FS) analysis === 
figure;

% Plot OS etc
PPOSsettings.PlotCPs = 0; % plot cardinal points (F,H,N)
PPOSsettings.PlotApertureAnalysis = [ 1 1 ]; % Pupils,Windows
PPOSsettings.CPs_Legend = 0; % Plot cardinal points (F,H,N) legend
PPOSsettings.ApA_Legend = [1 1]; % Pupils+Windows legend
PlotParaxialOS( OS , n12 , WD , PPOSsettings ); 

% Define a bunch of rays targetted to the axial-point of entrance pupil,
% with increasing angle (tilt) -- they will pass from AS axial point
Nr = 20;
rays = NaN*ones(Nr,4);
rays(:,1) = -WD; % S-position (axial)
rays(:,3) = 30 * linspace(-1,1,Nr) * pi/180; % angle [rad]
rays(:,2) = -rays(:,3).*(Pupils.enS+WD); % u-position (height)
wls = round( linspace( 633, 455, Nr/2) )'; % wavelength (color) [nm]
rays(:,4) = [wls; flipud(wls)];

% Do ray-tracing and plot rays
sua = DoParaxialRayTracing( OS, n12, rays, sfinal, 0 , 1);
hold on;
for kr = 1 : Nr
    plot( sua(:,1,kr), sua(:,2,kr), '-' , 'Color' , ...
        calculateVisibleSpectrumColor(rays(kr,4))  )
end

% Define and plot limiting rays, from object-plane axial-point to P_en edges
AFV = Windows.enD/(-Windows.enS+Pupils.enS);
raysFS = NaN*ones(2,4);
raysFS(:,1) = Windows.enS; % S-position (axial)
raysFS(:,2) = -Windows.enD*[-1/2,+1/2]; % u-position (height)
raysFS(:,3) = AFV * [-1/2 +1/2]; % angle [rad]
raysFS(:,4) = 700; % wavelength (color) [nm]
suaFS = DoParaxialRayTracing( OS, n12, raysFS, sfinal, 0 , 1);
for kr = 1 : size(raysFS,1)
    plot( suaFS(:,1,kr), suaFS(:,2,kr), '--' , 'Color' , ...
        calculateVisibleSpectrumColor(raysFS(kr,4))  );
end

title( 'Windows / Field-Stop analysis' );

