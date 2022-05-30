function PlotParaxialOS( OS , n12 , WD , PPOSsettings ) 

% Plots Paraxial OS (optical surfaces & materials) together with analysis
% annotations (cardinal points, apertures analysis etc).
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
    
    n12 = [1.0 1.3];
    WD = 200;
end


%% Set Options (flags,colors etc)

% Plot settings (defaults)
if nargin < 4 
    PPOSsettings.PlotCPs = 1; % plot cardinal points (F,H,N)
    PPOSsettings.PlotApertureAnalysis = 0*[ 1 1 ]; % Pupils,Windows
    PPOSsettings.CPs_Legend = 1; % Plot cardinal points (F,H,N) legend
    PPOSsettings.ApA_Legend = 0*[1 1]; % Pupils+Windows legend
end
PlotCPs = PPOSsettings.PlotCPs;
PlotApertureAnalysis = PPOSsettings.PlotApertureAnalysis; 
CPs_Legend = PPOSsettings.CPs_Legend; 
ApA_Legend = PPOSsettings.ApA_Legend; 

% Colors for filled pathces (lens body)
cglass = 0.85*[1 1 1];%  hex2dec( reshape('5580b6',[2,3])' )'/255; %
nglass = max( [OS(:,3);n12(:)] ); % reference index for color above
cair   = [1 1 1];%   hex2dec( reshape('4169aa',[2,3])' )'/255 ;% 

% Colors for cardinal points
cF = [0 0.5   1 ]; % "F" for First (front) focal, principal, nodal
cS = [1   0.5 0 ]; % "S" for Second (back) focal, principal, nodal

% Colors for lines
cP = [0 .5 0 ]; % Pupils
cW = [.5 0  .5]; % Windows


%% Misc/derived params

% Legend/display stuff:
LegHandles = [];
LegStrings = {};

% Retrieve paraxial analysis stuff
[~,CPs,Pupils,Windows] = DoParaxialAnalysis( OS , n12 , WD );

%Is afocal (F/H/N to infinity, or near-infinity), e.g. as in telescope?
EFL_OS = CPs.H1 - CPs.F1; % test EFL only
dOAT = OS(end,1) - OS(1,1); % physical length of OS
IsTelescope = all( abs( EFL_OS ) > 1000*dOAT ) && dOAT~=0; 

% Size of plot window (vertical)
if all( isinf( OS(:,4) ) )
    rmax = dOAT/2;
else
    auxi = OS(:,4);
    if PlotApertureAnalysis(1) == 1, auxi = [auxi;Pupils.enD]; end
    if PlotApertureAnalysis(2) == 1, auxi = [auxi;Windows.enD]; end
    rmax = 1.2*max( abs( auxi/2 ) ); 
end

% z-limits (used for plotting front/back materials, if n12>1)
z1 = OS(1,1) - WD; % limit in front
z2 = OS(end,1) + WD; % limit in back

% Assume all elements are lenses (or AS's if P=0)
if size(OS,2) == 4
    OS = [ OS , ones(size(OS,1),1) ]; 
end

% Prepare figure/axis
hold on;


%% Filled rectangles (lens bodies etc).

% Background color (for void, materials w/ refr.index == 1)
set(gca, 'Color' , cair )

% Plot "filled" glass areas (Rectangles)
ncurr = n12(1); % entrance medium
z1a = z1;
for kos = 0:size(OS,1)-1 % entrace material & all intermediate ones
    z2a = OS(kos+1,1);
    if ncurr ~= 1
        col = (cglass - cair)/(nglass - 1 + eps)*(ncurr - 1) + cair;
        fill( [z1a z1a z2a z2a] , rmax*[-1 +1 +1 -1] , col , ...
            'EdgeColor' , 'none' , 'FaceAlpha' , 1 )
    end
    z1a = OS(kos+1,1);
    ncurr = OS(kos+1,3);
end
if n12(2) ~= 1 % exit material
    col = (cglass - cair)/(nglass - 1 +eps)*n12(2) + cair;
    fill( [z1a z1a z2 z2] , rmax*[-1 +1 +1 -1] , col , ...
        'EdgeColor' , 'none' , 'FaceAlpha' , 1 )
end


%% Optical Surfaces (OS, verticals) ...
%La = 1; % length of arrow of lens (converging/diverging)
for kos = 1:size(OS,1)
    
    % ... depending on its Type:
    switch OS(kos,5)
        
        case 0, % Films (i.e. free-standing stops)
            plot( OS(kos,1)*[1 1], OS(kos,4)*[-1/2 +1/2], 'Color' , 'k' ,...
                'LineWidth' , 2 )
        
        case 1, % Lenses (w/ Stops) and Aperture Stops
            % Plot aperture-stop (up & down plate), for stops & lenses (assumed)
            if OS(kos,2)~=0
                col = 0.0*[1 1 1]; linwid = 1;
            else
                col = 'k'; linwid = 2;
            end
            plot( OS(kos,1)*[1 1], -[OS(kos,4)/2 rmax], 'Color' , col , 'LineWidth' , linwid )
            plot( OS(kos,1)*[1 1], +[OS(kos,4)/2 rmax], 'Color' , col , 'LineWidth' , linwid )
            
            % Check type (conv/div/stop): from value of optical power
            sgn = sign(OS(kos,2));
            if OS(kos,2)~=0 % refr surf: plot its surface and conv/div arrows
                col = [0 0 1];
                plot( OS(kos,1)*[1 1], OS(kos,4)/2*[-1 +1], 'Color' , col, 'LineWidth' , 2 )
                %plot( OS(kos,1)+[-La 0 La] , +OS(kos,4)/2+sgn*La*[-1 0 -1] , ...
                %    'Color', col, 'LineWidth' , 2 )
                %plot( OS(kos,1)+[-La 0 La] , -OS(kos,4)/2-sgn*La*[-1 0 -1] , ...
                %    'Color', col, 'LineWidth' , 2 )
                
                if sgn > 0, mu='^'; md='v'; else md='^'; mu='v';  end
                plot( OS(kos,1) , +OS(kos,4)/2 , 'Color', col, 'Marker', mu ,...
                    'MarkerSize' , 10 , 'LineWidth' , 2 , 'MarkerFaceColor' , col)
                plot( OS(kos,1) , -OS(kos,4)/2 , 'Color', col, 'Marker', md, ...
                    'MarkerSize' , 10 , 'LineWidth' , 2 , 'MarkerFaceColor' , col )
                
            else % stop: plot an "x" on the OA on the position of stop
                plot( OS(kos,1), 0 , 'k+' , 'MarkerSize' , 5 )
                % Stops: plot small "ears" (left-right)
%                 plot( OS(kos,1)+[-La 0 La] , +OS(kos,4)/2+[ 0 0 0 ] , ...
%                     'Color' ,col )
%                 plot( OS(kos,1)+[-La 0 La] , -OS(kos,4)/2-[ 0 0 0 ] , ...
%                     'Color' ,col )
            end
    
        case 2, % Reflector (free standing)
            col = 0.5*ones(1,3);
            plot( OS(kos,1)*[1 1], OS(kos,4)*[-1/2 +1/2], ...
                'Color' , col , 'LineWidth' , 2 );
            
            sgn = sign(OS(kos,2));
            plot( OS(kos,1)*[1 1], OS(kos,4)/2*[-1 +1], 'Color' , col, 'LineWidth' , 2 )
            %plot( OS(kos,1)+[-La 0 La] , +OS(kos,4)/2+sgn*La*[-1 0 -1] , ...
            %    'Color', col, 'LineWidth' , 2 )
            %plot( OS(kos,1)+[-La 0 La] , -OS(kos,4)/2-sgn*La*[-1 0 -1] , ...
            %    'Color', col, 'LineWidth' , 2 )
            
            if sgn > 0, mu='^'; md='v'; else md='^'; mu='v';  end
            plot( OS(kos,1) , +OS(kos,4)/2 , 'Color', col, 'Marker', mu ,...
                'MarkerSize' , 10 , 'LineWidth' , 2 , 'MarkerFaceColor' , col)
            plot( OS(kos,1) , -OS(kos,4)/2 , 'Color', col, 'Marker', md, ...
                'MarkerSize' , 10 , 'LineWidth' , 2 , 'MarkerFaceColor' , col )
                      
        case 3, % Reflector w/ aperture
            col = 0.5*ones(1,3);
            plot( OS(kos,1)*[1 1], -[OS(kos,4)/2 rmax], 'Color' , col , 'LineWidth' , 2 )
            plot( OS(kos,1)*[1 1], +[OS(kos,4)/2 rmax], 'Color' , col , 'LineWidth' , 2 )
            
           sgn = sign(OS(kos,2));
            plot( OS(kos,1)*[1 1], OS(kos,4)/2*[-1 +1], 'Color' , col, 'LineWidth' , 2 )
            %plot( OS(kos,1)+[-La 0 La] , +OS(kos,4)/2+sgn*La*[-1 0 -1] , ...
            %    'Color', col, 'LineWidth' , 2 )
            %plot( OS(kos,1)+[-La 0 La] , -OS(kos,4)/2-sgn*La*[-1 0 -1] , ...
            %    'Color', col, 'LineWidth' , 2 )
            
            if sgn > 0, mu='^'; md='v'; else md='^'; mu='v';  end
            plot( OS(kos,1) , +rmax , 'Color', col, 'Marker', mu ,...
                'MarkerSize' , 10 , 'LineWidth' , 2 , 'MarkerFaceColor' , col)
            plot( OS(kos,1) , -rmax , 'Color', col, 'Marker', md, ...
                'MarkerSize' , 10 , 'LineWidth' , 2 , 'MarkerFaceColor' , col )
            
    end
    
end


%% Cardinal points (when P~=0)
if ~IsTelescope && PlotCPs
    cph = NaN*ones(1,6);
    cph(1) = plot( CPs.F1*[1 1] , 0*[-1 1] , 'o' , 'Color' , cF, 'MarkerFaceColor' , cF , 'MarkerSize' , 6 );
    cph(2) = plot( CPs.F2*[1 1] , 0*[-1 1] , 'o' , 'Color' , cS, 'MarkerFaceColor' , cS , 'MarkerSize' , 6 );
    cph(3) = plot( CPs.H1*[1 1] , rmax*[-1 1], ':' , 'Color' , cF , 'LineWidth' , 1);
    cph(4) = plot( CPs.H2*[1 1] , rmax*[-1 1], ':' , 'Color' , cS , 'LineWidth' , 1);
    cph(5) = plot( CPs.N1*[1 1] , 0*[-1 1] , '+' , 'Color' , cF , 'MarkerSize' , 12 );
    cph(6) = plot( CPs.N2*[1 1] , 0*[-1 1] , '+' , 'Color' , cS , 'MarkerSize' , 12 );
    
    %set(gca,'Layer','bottom')
    if CPs_Legend == 1
        LegHandles = [ LegHandles , cph ];
        LegStrings = [ LegStrings , {'F_1','F_2','H_1','H_2','N_1','N_2'} ];
    end
end



%% Apertures (AS/FS, Pupils & Windows) 

% Plot Entrance Pupil (and mark AS)
if PlotApertureAnalysis(1) == 1
    if ~isempty( Pupils.iAS ) && abs(Pupils.enS) < 1000*dOAT 
        
        if Pupils.enS > OS(1,1)
            penls = '--'; % P_entrance_line_style
        else
            penls = '-';
        end
        
        hps = NaN*ones(1,3);
        hps(1) = plot( Pupils.enS*[1 1], +1*[abs(Pupils.enD)/2 rmax] , penls ,...
            'Color' , cP, 'Linewidth' , 2 );
        plot( Pupils.enS*[1 1], -1*[abs(Pupils.enD)/2 rmax] , penls , ...
            'Color' , cP, 'Linewidth' , 2 )

        plot( Pupils.enS , 0 , '+' , 'Color' , cP , ...
            'MarkerFaceColor' , cP , 'MarkerSize' , 8 )

        hps(2:3) = plot( OS(Pupils.iAS,1) , OS(Pupils.iAS,4)*[-1/2 +1/2] , 'o' , ...
            'Color' , cP , 'MarkerFaceColor' , cP , 'MarkerSize' , 8 );
        
        if ApA_Legend(1) == 1
            LegHandles = [ LegHandles , hps(1:2) ];
            LegStrings = [ LegStrings , {'P_{en}','AS'} ];
        end
                
    end
end

% Plot Entrance Window (and mark FS)
if PlotApertureAnalysis(2) == 1
    if ~isempty( Windows.iFS ) && abs(Windows.enS) < 1000*dOAT
        
        if Windows.enS > OS(1,1)
            wenls = '--'; % W_entrance_line_style
        else
            wenls = '-';
        end
        
        hws = NaN*ones(1,3);
        hws(1) = plot( Windows.enS*[1 1], +1*[abs(Windows.enD)/2 rmax] , wenls , ...
            'Color' , cW, 'Linewidth' , 4 );
        plot( Windows.enS*[1 1], -1*[abs(Windows.enD)/2 rmax] , wenls, ...
            'Color' , cW, 'Linewidth' , 4 )

        plot( Windows.enS , 0 , '+' , 'Color' , cW , ...
            'MarkerFaceColor' , cW , 'MarkerSize' , 8 )

        hws(2:3) = plot( OS(Windows.iFS,1) , OS(Windows.iFS,4)*[-1/2 +1/2] , 'o' , ...
            'Color' , cW , 'MarkerFaceColor' , cW , 'MarkerSize' , 8 );
        
        if ApA_Legend(2) == 1
            LegHandles = [ LegHandles , hws(1:2) ];
            LegStrings = [ LegStrings , {'W_{en}','FS'} ];
        end
        
    end
end

legend( LegHandles , LegStrings )


%% Optical Axis (OA) & Scaling
zLim = get(gca,'XLim') + [ -10 +10 ];
zLim(1) = min( [zLim(1), -WD] );
plot( zLim , [0 0], 'm' , 'LineWidth', 1 ) % optical axis (OA)
plot( max(zLim) , 0 , 'm>', 'MarkerFaceColor' , 'm' , 'MarkerSize', 10);
set( gca, 'YLim' , rmax*[-1 +1] )

set(gca,'Layer','top')
xlabel( 'S-axis (absolute)')
ylabel( 'u-axis')

% Toggle-button for axis scaling
axis equal;
htas=uicontrol('Style','togglebutton','Position',[10 10 80 25],...
    'String' , 'Scale (equal)' , 'Value' , 1 , 'Callback', @toggleAxisEqual );
set(gcf,'WindowStyle','normal')
set(gcf,'Toolbar','figure')

end

% GUI Function: Toggle axis scaling (equal/off)
function toggleAxisEqual(htas,eventdata)
    switch get(htas,'Value')
        case 1, 
            %axis equal;
            set(gca,'dataaspectratio',[1 1 1]); 
            set(htas,'String','Scale (equal)')
        case 0, 
            axis normal; 
            %Lx = abs(diff( get(gca,'XLim') ) );
            %Ly = abs(diff( get(gca,'YLim') ) );
            %set(gca,'dataaspectratio',[Ly/Lx, 1, 1] );
            set(htas,'String','Scale (off)')
    end
end


