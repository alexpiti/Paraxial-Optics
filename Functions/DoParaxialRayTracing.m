function sua = DoParaxialRayTracing( OS, n12, rays, ...
    sfinal, doPlot , EnableApertureStops )

% Paraxial ray tracing, given an optical system defined with the OS
% (optical surfaces) matrix, the n12 (entrace/exit refr.indx) and a set of
% rays.
%
% *** Lengths are normalized. Just be careful to use consistent units
%     between lengths and optical powers (1/length)
% 
% >> OS (Optical Surfaces)
%    A Nos-by-4 (or Nos-by-5) array containing the information of Nos
%    optical surfaces (paraxial==plane). Parameters in the columns of OS as:
%  - 1. Absolute s-position (z-horizontal), with surfaces & stops 
%       assumed seq. arranged. OA folding allowed.
%  - 2. Optical power (units: Length^-1). P=0 denotes an Aperture Stop (AS)
%  - 3. Refractive index to material on its RIGHT (exit)
%  - 4. Diameter (assumed cylindrical-symmetric & centered on OA). This is
%       required for apertures (pupils & windows) analysis.
%  * 5. [OPTIONAL] Type of optical surface: 0:stop (free-standing absorber)
%       {1:Lens or AS}, 2:mirror(free-standing), 3:mirror (with an aperture
%       in its middle). In this implementation, Type==1 is the default. 
%       The others are needed for modeling catadioptric structures.
%
% >> rays
%    A Nr-by-4 array containing the information of the total of Nr rays. 
%    The parameters of the ray are in its four columns as:
%  - 1. Absolute S-position of ray source
%  - 2. u-position of ray source
%  - 3. angle of ray (radians)
%  - 4. wavelength of ray (in nm), used for its color
%
% >> sua
%  Stands for "s-u-a" = {S-position, u-position, angle} of a ray. This is a
%  (Nos+2)-by-3-by-Nr array, containing all the information of the rays
%  passage through the optical system.
% 
% GNU General Public License v3.0 | Copyright (c) 2022 Alexandros Pitilakis
% (c) Alexandros Pitilakis, Thessaloniki/Greece, April 2019


%% Test inputs
if nargin == 0
    clc; close all; hold on;
    
    OS(1,:) = [ 0   1/20 1.50 20 1 ];
    OS(2,:) = [ 5      0 1.50 10 1 ]; % acting as AS
    OS(3,:) = [ 10 -1/20 1.00 20 1 ];
    OS(4,:) = [ 70     0 NaN  30 1 ]; % acting as FS (e.g. film)
    
    n12 = [1 2];
    WD  = 200;
    
    PlotParaxialOS( OS , n12 , WD );
        
    % Define a few rays
    S1 = -WD;
    rays(1,:) = [ S1, 5, +0.07, 633];
    rays(2,:) = [ S1, 5,  0,    550];
    rays(3,:) = [ S1, 5, -0.07, 455];
    
    % Optional inputs
    sfinal = 2*OS(end,1);
    doPlot = 1;
    EnableApertureStops = true;
    
end


%% Function Flags & Parameters
DrawFilledRayCone = 0;


%% Prepare matrices and variables

if nargin < 4, sfinal = []; end
if nargin < 5, doPlot = true; end
if nargin < 6, EnableApertureStops = false; end    

% Assume all elements are lenses (AS)
if size(OS,2) == 4
    OS = [ OS , ones(size(OS,1),1) ]; 
end

% Retrieve ray parameters
Nr = size(rays,1); % number or rays
sr = rays(:,1); % s-positions (starting) of rays
ur = rays(:,2); % u-positions (staring) of rays
ar = rays(:,3); % angles (starting) of rays
cr = rays(:,4); % colors of rays (wavelength in nm)

% OS and stuff
[ABCD] = DoParaxialAnalysis( OS , n12 );

% All-ray data, i.e s-u-a (s:horiz, u:vert a:angle[rad])
% ** Nrows = 1(input) + N_opticalSurfaces + 1(output) = N+2
sua = NaN*ones(size(OS,1)+2,3,Nr);


%% Sweep rays in bundle
for kr = 1:Nr
    
    % S-final
    if isempty(sfinal)
        sfinal = DoParaxialImageFormation( ABCD , n12 , OS(end,1)-OS(1,1) , ...
            OS(1,1) , sr(kr) , [] ) % absolute s-pos of image
        sfinal = max( [sfinal,get(gca,'XLim')] );  
    end
    
    % Initializations
    stopped = false; % does this ray terminate on a stop?
    n = n12(1); % current material refr.indx == entrance material
    clear una
    sua(1,:,kr) = [sr(kr), ur(kr), ar(kr)]; % start pos (su) of this ray
    una = [ ur(kr) ; n*ar(kr) ]; % [ (vert.pos) ; (refr.indx) * (angle_rad) ]
    
    % Original direction of rays (assumed to OA direction, +). This is
    % required for reflective OS, which change the direction of
    % S-propagation.
    Ssgn = +1;    
    
    % Sweep all optical surfaces (assumed sequentially arranged)
    % ABCD matrices are assumed
    for kos = 1:size(OS,1)
        
        % Propagate up to optical surface
        una = [1 (OS(kos,1)-sua(kos,1,kr))/n*Ssgn; 0 1]*una;
        
        % Store ray position & angle (on optical surface)
        sua(kos+1,:,kr) = [OS(kos,1), una(1), una(2)/n];
        
        % Check whether ray is blocked by an aperture stop
        if EnableApertureStops == false || ...
          ((OS(kos,5)==0 && abs(una(1)) >  OS(kos,4)/2) || ... % Stop (free-standing)
           (OS(kos,5)==1 && abs(una(1)) <= OS(kos,4)/2) || ... % Lens or Ape
           (OS(kos,5)==2 && abs(una(1)) <= OS(kos,4)/2) || ... % Mirror (free-stand)
           (OS(kos,5)==3 && abs(una(1)) >  OS(kos,4)/2))        % Mirror w/ aperture
            
           if OS(kos,5)>=2 % reflective
               Ssgn = -1*Ssgn;
           end
       
            % Refract on OS (and proceed)
            una = [1 0; -OS(kos,2) 1]*una;
            
            % Store ri of material after (right) of opt-surf
            n = OS(kos,3);
            
        else
            stopped = true;
            break
        end
        
    end
    
    % if ray wasn't stopped, it can proceed until the assumed end of OS
    if stopped == false
        % One last propagation up to the final position
        una = [1 (sfinal-OS(kos,1))/n12(2);0 1]*una;
        
        % Store last ray position & angle
        sua(kos+2,:,kr) = [sfinal, una(1), una(2)/n12(2)];
    end
    
    % Plot this ray from is "s-u-a" data
    if doPlot == 1
        plot( sua(:,1,kr) , sua(:,2,kr) , 'Color' , ...
            calculateVisibleSpectrumColor(round(cr(kr))) , ...
            'Linewidth' , 2 )    
    end
    
end

% Ray cone
if doPlot == true && Nr > 1 && DrawFilledRayCone == true
    i1 = find( squeeze(sua(1,3,:)) == max( sua(1,3,:) ) );
    i2 = find( squeeze(sua(1,3,:)) == min( sua(1,3,:) ) );
    xs = [sua(:,1,i1); flipud(sua(:,1,i2)) ];
    ys = [sua(:,2,i1); flipud(sua(:,2,i2)) ];
    patch( xs, ys, calculateVisibleSpectrumColor(round(cr(kr))) , ...
       'EdgeColor' , 'None' , 'FaceAlpha' , 0.3 )
end

