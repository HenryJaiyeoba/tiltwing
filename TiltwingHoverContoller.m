function phiwingset = TiltwingHoverContoller(theta_err,  dthetadt_err, phiwing0, h, l)

%PD-controller: alpha = Kd_tild * dthetadt_err + Kp_tild * theta_err,
%where alpha - tiltwing deflection from neutral orientation, theta_err - UAV
%pitch error

%INPUTS
%phiwing0 - at this angle the moment induced by prop thrust must be zero (or at least close to zero)
%h - leg between CG and tiltwing axis
%l - characteristic length of UAV (inertia moment I = m * l^2)

g = 9.81;

Kp_tild = 0.8; %how many degrees is phi when theta is 1 degree - this is the main tunable parameter

Kp = g*h/l^2 * Kp_tild;
Kd = 2*sqrt(Kp);
Kd_tild = Kd / (g*h/l^2);

alpha = Kd_tild * dthetadt_err + Kp_tild * theta_err;
phiwingset = phiwing0 - alpha;
