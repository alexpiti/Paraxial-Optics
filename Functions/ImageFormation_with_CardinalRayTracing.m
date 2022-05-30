function [ABCD,CPs,imaS,imau] = ImageFormation_with_CardinalRayTracing(...
    OS, n12, objS, obju, sfinal )

% Draw cardinal rays from given vertical object (objS,obju), real or
% virtual, through the given OS. n12 are the indices of materials to the
% entrace and exit of the OS.
%
% *** objS is the ABSOLUTE S-position of the object. Typically, negative
%     S-value denotes a real object (before the entrance).
%
% GNU General Public License v3.0 | Copyright (c) 2022 Alexandros Pitilakis
% (c) Alexandros Pitilakis, Thessaloniki/Greece, April 2019

% Perform paraxial analysis (for cardinal points, CPs)
[ABCD,CPs] = DoParaxialAnalysis( OS , n12 );

% Analytically calculate the (Guassian) image position & magnification
[imaS,~,mu] = DoParaxialImageFormation( ABCD , n12 , OS(end,1)-OS(1,1) ,...
    OS(1,1) , objS , [] );
imau = mu*obju;


% Plot the OS and annotations (cardinal points etc)
PlotParaxialOS( OS , n12 , abs(objS) ); 

% Rays will be propagated up to this s-position
if nargin < 5 || isempty( sfinal )
    sfinal = 2*max( abs([imaS,objS,OS(end,1)]) );
end

% Rays are defined in the "rays" matrix. Each row is a different ray and
% the parameters of each ray are:
%   column #1: Absolute S-position of ray source
%   column #2: Transverse (u) position of ray source
%   column #3: Angle of ray wrt OA (radians)
%   column #4: wavelength (nm) of ray, for coloring it.
%
% In image formation, we define three "cardinal" rays:
if objS < OS(1,1) % Real object (i.e. all rays START from object)
    rays(1,:) = [ objS, obju, 0, 633]; % parallel to OA in entrance
    rays(2,:) = [ objS, obju, (0-obju)/(CPs.F1-objS), 633 ]; % pointing to F1
    rays(3,:) = [ objS, obju, (0-obju)/(CPs.N1-objS), 633 ]; % pointing to N1
else % virtual object (rays DO NOT cross in real object space)
    comS = min( [-2*abs(objS),-sfinal] ); % common S-start for all rays
    dS = comS - objS; % S-distance of comS with virtual object 
    rays(1,:) = [ comS, obju, 0, 633]; % parallel to OA
    angF1 = obju/(objS-CPs.F1); uF1 = obju + angF1*dS;
    rays(2,:) = [ comS, uF1 , angF1 , 633 ];
    angN1 = obju/(objS-CPs.N1); uN1 = obju + angN1*dS;
    rays(3,:) = [ comS, uN1 , angN1 , 633 ];
end

% Do ray-tracing (plots rays)
sua = DoParaxialRayTracing( OS, n12 , rays , sfinal );

% ------------------------------------------

% Plot object (real or virtual)
if objS > OS(1,1), linst = '--'; else linst = '-'; end
col = [1 0 0];
ho=plot(  objS*[1 1] , obju*[0 1] , 'Linewidth' , 4 , 'Color' , col , 'LineStyle' , linst );
if obju > 0, marker = '^'; else marker = 'v'; end
plot(  objS, obju , 'Marker' , marker , 'MarkerSize', 10, ...
    'MarkerFaceColor' , col , 'Color' , col   );

% Trace-back virtual rays from virtual object to entrance surface (A1)
if objS > OS(1,1)
    for kr = 1:3
        plot( [objS,sua(2,1,kr)] , [obju,sua(2,2,kr)], 'r:' )
    end
end
% Trace virtual rays from object to F1 & N1
plot( [objS,CPs.F1] , [obju,0], 'r:' )
plot( [objS,CPs.N1] , [obju,0], 'r:' )


% Plot the analytically calculated image (real or virtual)
if imaS < OS(end,1), linst = '--'; else linst = '-'; end
col = [0 0.5 0];
hi = plot(  imaS*[1 1] , imau*[0 1] , 'Linewidth' , 4 , 'Color' , col , 'LineStyle' , linst );
if imau > 0, marker = '^'; else marker = 'v'; end
plot(  imaS, imau , 'Marker' , marker , 'MarkerSize', 10, ...
    'MarkerFaceColor' , col , 'Color' , col  );

% Trace-ahead rays from virtual image to exit surface (A2)
if imaS < OS(end,1)
    for kr = 1:3
        plot( [imaS,sua(end-1,1,kr)] , [imau,sua(end-1,2,kr)], ':',...
            'color', [0 .5 0])
    end
end
% Trace virtual rays from image to F2 & N2
%plot( [imaS,CPs.F2] , [imau,0], ':' , 'Color' , col )
%plot( [imaS,CPs.N2] , [imau,0], ':' , 'Color' , col )

% Re-plot OA
delete( findobj('Color','m') )
plot( get(gca,'XLim')+10*[-1 +1] , [ 0 0 ] , 'm' )
plot( max(get(gca,'XLim')) , 0 , 'm>', 'MarkerFaceColor' , 'm' , 'MarkerSize', 10);

% Grab current legend user-data (for existing hangles & strings)
try
    clud = get( legend, 'UserData' );
    clh = clud.handles;
    cls = clud.lstrings;
    legend( [clh;ho;hi], [cls;'Object';'Image'] )
catch
    legend( [ho;hi], {'Object';'Image'} )
end