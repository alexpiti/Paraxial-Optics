function OS = Convert_ABCD2OS( ABCD , d , S0 )

% Converts ABCD matrix to equivalent thick (or thin) lens, i.e., in an
% Optical Surface (OS) representation. S0 is the absolute S-position of the
% entrance vertex (A1) of the system.
%
% GNU General Public License v3.0 | Copyright (c) 2022 Alexandros Pitilakis
% (c) Alexandros Pitilakis, Thessaloniki/Greece, April 2019

% Checked determinant of ABCD.
if abs(det(ABCD)-1)>10*eps
   error( ' Invalid ABCD: Check that det(ABCD)==1' );
end

A = ABCD(1);
C = ABCD(2);
B = ABCD(3);
D = ABCD(4);

% Represent with a thick (or thin) lens of thickness=d and appropriate 
% refractive index and optical powers (e.g. curvatures) in each side.
if d > 0 % Thick lens
    n  = d/B;
    P1 = (1-A)*n/d;
    P2 = (1-D)*n/d;
    OS(1,:) = [ S0   , P1 , n   , 5*d ];
    OS(2,:) = [ S0+d , P2 , NaN , 5*d ];
else 
    OS(1,:) = [ S0   , -C , NaN , abs(1/C) ] ; % Thin lens.
end

