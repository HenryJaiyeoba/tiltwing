function [omega_servo, alpha_servo] = ConfineOmegaServo(omega_set, omega_prev, dt, domega_dt_UAVz, omegamax, Mservomax, Jduct3)
%this function imitates natural limitations of servomotor, confining
%instant chance in angular speed
%! here omega is signed: positive - "duct up"

%SOME INPUTS:
%domega_dt_UAVz - rotation acceleration of UAV about OZ (nose up if positive)
%Mservomax - maximal torque of the duct servo
%omegamax - maximal servo rot. speed
%Jduct3 - inertia moment of the whole duct about OZ

d_omega = omega_set - omega_prev;
Mservomax_posit = Mservomax - Jduct3*domega_dt_UAVz; %which "duct up" motion servo can support
Mservomax_negat = Mservomax + Jduct3*domega_dt_UAVz; %which "duct down" motion servo can support

Mreq = Jduct3 * d_omega/dt;

Mreq(Mreq == 0) = 1e-9; %protection against /0

omega_servo = omega_prev + d_omega .* (min(1, Mservomax_posit./abs(Mreq)) .* (Mreq>0) + min(1, Mservomax_negat./abs(Mreq)) .* (Mreq<0));

omega_servo = max(-omegamax, min(omegamax, omega_servo));
alpha_servo = (omega_servo - omega_prev)/dt;