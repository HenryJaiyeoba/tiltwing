%function debug and tuning
close all

alpvec = pi * (-1:0.001:1); %wing angle of attack

Clvec = NaN(size(alpvec));
Cdvec = NaN(size(alpvec));
Fl_arr = NaN(2, length(alpvec));
Marr = NaN(size(alpvec));
Fn_arr = NaN(2, length(alpvec));

Rl_ac = [0; 0.2];
rho = 1.22;
S = 1.5;
chord = 0.4;

Clmin = -1.2;
Clmax = 1.4;
Cdmin = 0.02;
Cdmax = 1.1;
Cm = -0.03;
alpmax = 15*pi/180;

for n = 1:length(alpvec)
    Vl_air = -[cos(-alpvec(n)); sin(-alpvec(n))];
    [Fl_arr(:,n), Marr(n), Clvec(n), Cdvec(n)] = FaeroMaero(Vl_air, Rl_ac, rho, S, chord, Cm, Clmin, Clmax, Cdmin, Cdmax, alpmax);
    Cnl = OrientMatrix2D(alpvec(n));
    Fn_arr(:,n) = Cnl * Fl_arr(:,n);
end



figure
subplot(311)
plot(180/pi*alpvec, Clvec); grid on
subplot(312)
plot(180/pi*alpvec, Cdvec); grid on
subplot(313)
plot(180/pi*alpvec, Clvec./Cdvec); grid on

figure
subplot(311)
plot(180/pi*alpvec, Fn_arr(1,:), 'b', 180/pi*alpvec, Fn_arr(2,:), 'r'); grid on
subplot(312)
plot(180/pi*alpvec, 180/pi * (atan2(Fn_arr(2,:), Fn_arr(1,:)) - alpvec)); grid on
subplot(313)
plot(180/pi*alpvec, Marr); grid on

%check the direction of aero. force
figure
for n = 1:length(alpvec)
    plot([0 cos(alpvec(n))], [0 sin(alpvec(n))]); grid on; hold on
    plot([0 Fn_arr(1,n)], [0 Fn_arr(2,n)]); hold off
    axis(2*[-1 1 -1 1])
    drawnow
end










