function [Fb, Maero] = FaeroMaero2D(Vb_air, rho, S, chord, Cm, Clmin, Clmax, Cdmin, Cdmax, alpmax)
%this function calculates aerodynamic force (2D vector) and torque (scalar,
%anticlock.) of a wing/stabilitor/aileron/fuselage in "local" (airfoil-attached)
%coordinate system; lift and drag coefs - just for debugging

%Vb_air - airspeed in RF "b"
%alpmax - angle of attack where maximal lift coef. is observed

alp = -atan(Vb_air(2) / Vb_air(1)); %wing angle of attack
Vair = norm(Vb_air);
if Vair < 1e-6 %protection against NaN
    alp = 0;
end

Cl_amp = (Clmax - Clmin)/2;
Cl_shift = (Clmax + Clmin)/2;
Cl_1 = Cl_shift + Cl_amp * sin(pi/2* alp/alpmax);
Cl_2 = (Clmax*(alp>0) + Clmin*(alp<0)) * max(1 - 0.5*((abs(alp) - alpmax)/alpmax)^2, 0);
Cd_amp = (Cdmax - Cdmin)/2;
Cd_shift = (Cdmax + Cdmin)/2;

Cl = Cl_1 * (abs(alp) < alpmax) + Cl_2 * (abs(alp) >= alpmax);

Cd = Cd_shift - Cd_amp*cos(2*alp);

aero_coef = rho/2 * S * Vair^2;

Fl = Cl*aero_coef;
Fd = Cd*aero_coef;
Maero = Cm*chord*aero_coef;

alp_2pi = -atan2(-Vb_air(2), -Vb_air(1)); %wing angle of attack (-pi,pi)
if Vair < 1e-6 %protection against NaN
    alp_2pi = 0;
end
Cab = OrientMatrix2D(alp_2pi); %trans. matrix from "local" to "air" (where flight is "horizontal")
Fb = Cab' * [-Fd; Fl];

