function [ABCD,CPs,Pupils,Windows] = DoParaxialAnalysis( OS , n12 , WD )

% Paraxial analysis of a system of dioptric optical surfaces (OS) 
% with limiting apertures. Calculates total transfer matrix (ABCD) and its
% cardinal points, and also Pupils & Windows.
%
% --- Inputs ---
%  > OS (Optical Surfaces): A N-by-4 matrix where first row is the ABSOLUTE
%    S-position of OS-vertex (units=L), second row is optical power 
%    (units 1/L), third row is refractive index to the RIGHT (exit) of OS 
%    and fourth row is the diameter of the lens aperture (units=L). The
%    number of optical surfaces is N. OS are assumed sequentially arranged.
%      ** Power=0 denotes a simple aperture stop.
%  > n1: The refractive index to the LEFT (entrace) of the first OS
%  > WD (Working Distance): Distance of object-plane from the first OS.
%    Using the signed S-convention, i.e., a positive WD is to the LEFT 
%    (entrance) of the first OS, which is the typical case for imaging 
%    systems. If the first OS vertex is at S_abs=0, then WD_Sabs = -WD.
%
% --- Outputs ---
%  > ABCD of total optical system
%  > Cardinal Points (CPs): F1,F2,H1,H2,N1,N2 -- in absolute S-pos
%  > Pupils : iAS (index in rows of OS), enS (abs), enD, exS (abs), exD
%  > Windows: iFS (index in rows of OS), enS (abs), enD, exS (abs), exD
%
% === Notes ===
% * Treatment of lengths (and relevant dimensions) is normalized, meaning
%   that it's inferred that you must define all lengths (e.g. thicknesses,
%   diameters, distances) in the same units (e.g. mm or inches) and optical
%   powers correspondingly (e.g. in 1/mm or 1/inches).
% * Angles are treated (internally) in radians. 
%
% GNU General Public License v3.0 | Copyright (c) 2022 Alexandros Pitilakis
% (c) Alexandros Pitilakis, Thessaloniki/Greece, April 2019

% Prelims
n1 = n12(1);
n2 = n12(2); 

% Calculate ABCD matrix of ENTIRE OS. In parallel, find entrace pupil
ABCD = eye(2); % initialize ABCD
SPens = NaN*ones(1,size(OS,1)); SPens(1) = 0;
DPens = NaN*ones(1,size(OS,1)); DPens(1) = OS(1,4);
for kos = 1 : size(OS,1)-1
    
    % (1) Refract on current OS and (2) propagate to next one.
    OptSurfDist = abs( OS(kos+1,1)-OS(kos,1) ); % abs --> For REFLECTING/FOLDED!
    ABCD = [1 OptSurfDist/OS(kos,3); 0 1]*[1 0; -OS(kos,2) 1]*ABCD ;
    
    % Find images (S-pos & Diam) of each aperture, to the left, assuming
    % that each next aperture (at this.S2=0) is the IMAGE, so that its
    % conjugate OBJECT is the candidate pupil (we need its S1 & 1/mu) 
    A = ABCD(1); B = ABCD(3);
    SPens(kos+1) = -(-B*n1/A) + OS(1,1); % Absolute S-position (from S1)
    DPens(kos+1) = 1/A*OS(kos+1,4); % Diameter of pupil 
    
end
ABCD = [1 0; -OS(end,2) 1]*ABCD; % last surface refraction/reflection

% Final compacted optical system parameters
A = ABCD(1); B = ABCD(3);
C = ABCD(2); D = ABCD(4);

% Cardinal points (absolute S-positions, from A1 and A2). Store in a struct
CPs.F1 = OS( 1 ,1) -( -D/C*n1 ); 
CPs.F2 = OS(end,1) +( -A/C*n2 );
CPs.H1 = OS( 1 ,1) -( n1*(1-D)/C ); 
CPs.H2 = OS(end,1) +( n2*(1-A)/C );
CPs.N1 = OS( 1 ,1) -( (n2-n1*D)/C ); 
CPs.N2 = OS(end,1) +( (n1-n2*A)/C );

% In case there's no limiting apertures...
if nargin < 3 || all( isinf( OS(:,4) ) )
    Pupils.iAS = [];
    Pupils.enS = [];
    Pupils.enD = [];
    Pupils.exS = [];
    Pupils.exD = [];
    Windows.iFS = [];
    Windows.enS = [];
    Windows.enD = [];
    Windows.exS = [];
    Windows.exD = [];
    return
end

% Locate limiting Apperture Stop and Entrance Pupil (P_en): Entrance
% pupil is the one "seen" w/ the least angle from the on-axis WD (working
% distance) point. Equivalently, it's the aperture that is the first to
% stop the oblique ray from the on-OA object (at WD).
if ~isinf(WD)
    AngCandPens = abs(0.5*DPens./(SPens+WD));
    iAS = find( AngCandPens == min(AngCandPens) , 1 ); % index of AS (in OS)
else
    iAS = find( DPens == min(DPens) , 1 );
end
SPen = SPens( iAS ); % Absolute S-position of THE entrance pupil
DPen = DPens( iAS ); % diameter

% Calculate exit pupil
S1 = OS(1,1) - SPen; % relative to A1 (physical entrace)
S2 = -n2*(B*n1+A*S1)/(D*n1+C*S1); % relative S-pos
SPex = S2 + OS(end,1); % absolute S-pos
DPex = DPen * ( A + C*S2/n2 ); % diameter

% Store in a structure
Pupils.iAS = iAS;
Pupils.enS = SPen;
Pupils.enD = DPen;
Pupils.exS = SPex;
Pupils.exD = DPex;

% Locate limiting Field Stop (FS) and Entrance Windox (W_en). Here we use
% the alternative method: Launch a ray targeting the P_en, with a very
% small tilt so that it will surely traverse all the OS (along the OA). The
% FS is the aperture where |u|/Diameter is maximal, with u being the 
% transverse distance of the ray to the OA.
if size(OS,1) >= 2
    
    % Test paraxial ray (ang<<1) directed to Pen center
    suw = NaN*ones( size(OS,1)+1 , 2); % each row: [s,u] of this "window" ray
    n = n1; % current material refr.indx.
    swr = OS(1,1)-10; % starting S-position (absolute) of this ray
    uwr = 0.01; % [rad] starting u-position of this ray
    awr = uwr / -( SPen - swr ); % target ray to entrance pupil
    suw(1,:) = [swr uwr]; % start pos (S,u) of this ray
    una = [ uwr ; n*awr ]; % [ (vert.pos) ; (refr.indx) * (angle_rad) ]
        
    % Sweep all optical surfaces (assumed sequentially arranged) and
    % calcuate u-pos on each one of them.
    for kos = 1:size(OS,1)        
        % Propagate up to optical surface
        una = [1 (OS(kos,1)-suw(kos,1))/n; 0 1]*una;        
        % Store ray position (on optical surface)
        suw(kos+1,:) = [OS(kos,1), una(1)];        
        % Refract on OS (and proceed)
        una = [1 0; -OS(kos,2) 1]*una;        
        % Store ri of material AFTER (right) of this OS
        n = OS(kos,3); 
    end
    
    % Calculate rations (|u|/Diameter) and find the maximal
    DiamRatios = abs(suw(2:end,2)) ./ OS(:,4);
    iFS = find( DiamRatios == max(DiamRatios) , 1 );
    
    % Calc Entrance Window (W_en) S-pos & Diameter
    if iFS == 1 % foremost OS-aperture is the FS ==> FS == Wen
        SWen = OS(1,1); % Absolute S-pos of Entrance Window
        DWen = OS(1,4); % Diameter of Wen
    else        
        % We need the ABCD matrix (here called "T") from the first OS up to
        % the OS right before the FS
        T = eye(2);
        for kos = 1:iFS-1
            % (1) Refract on OS, (2) propagate
            OptSurfDist = abs( OS(kos+1,1)-OS(kos,1) ); % abs --> For REFLECTING/FOLDED!
            T = [1 OptSurfDist/OS(kos,3); 0 1]*[1 0; -OS(kos,2) 1]*T ;
        end
        SWen = -(-T(3)*n1/T(1)) + OS(1,1); % Absolute W_ex S-pos (from S1, for S2=0)
        DWen = 1/T(1)*OS(iFS,4); % Diameter of Wen   
    end
    
    % Calculate exit pupil
    S1 = OS(1,1) - SWen; % relative to A1 (physical entrace)
    S2 = -n2*(B*n1+A*S1)/(D*n1+C*S1); % relative S-pos
    SWex = S2 + OS(end,1); % absolute S-pos
    DWex = DWen * ( A + C*S2/n2 ); % diameter
    
else % For a single OS, there is no FS (AFV==180deg)
    
    iFS  = [];
    SWen = [];
    DWen = [];    
    SWex = [];
    DWex = [];    

end

% Store in a structure
Windows.iFS = iFS;
Windows.enS = SWen;
Windows.enD = DWen;
Windows.exS = SWex;
Windows.exD = DWex;
