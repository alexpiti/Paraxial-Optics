function [S2,S1,mu] = DoParaxialImageFormation( ABCD , n12 , d , S0 , S1 , S2 )

% Paraxial image formation.
%
% Given the ABCD & n12 (S0 is the absolute S-position of the A1: front
% physical vertex of the OS; d is the thickness of the OS, so that the back
% physical vertex is A2 at S0+d), and ONE of S1 or S2 (absolute positions), 
% this function returns the S2 or S1 (absolute positions) and the 
% magnification.
%
% GNU General Public License v3.0 | Copyright (c) 2022 Alexandros Pitilakis
% (c) Alexandros Pitilakis, Thessaloniki/Greece, April 2019

% Test inputs
if nargin == 0
   ABCD =  [0.6 3.8; -0.2 0.4]; 
   n12 = [1 2]; d = 10; S0 = 0;
   S1 = [];
   S2 = -34;
end

% Error check:
if abs(det(ABCD)-1) > 100*eps,
    error( ' ABCD matrix must satisfy: det(ABCD) == 1' );
end

% Grab parameters
A = ABCD(1); B = ABCD(3);
C = ABCD(2); D = ABCD(4);
n1 = n12(1); n2 = n12(2);

% Check if OBJECT or IMAGE S-pos is given
if ~isempty(S1) && isempty(S2) % Calc image position
    S1r = S0 - S1; % relative to A1 = physical entrance
    S2r = -n2*(B*n1 + A*S1r)/(D*n1 + C*S1r); % relative to A2 = phys. exit
    S2  = S0 + d + S2r; % absolute   
    mu  = A + C*S2r/n2; % magnification
elseif ~isempty(S2) && isempty(S1)
    S2r = S2 - S0 - d; % relative to A2
    S1r = -n1*(B*n2 + D*S2r)/(A*n2 + C*S2r); % relative to A1
    S1  = S0 - S1r; % absolute
    mu  = 1/(A+C*S2r/n2); % magnification
else
    error(' Please provide only S1 or only S2' )
end

% Test outputs
if nargin == 0
   S1
   S2
   mu
end