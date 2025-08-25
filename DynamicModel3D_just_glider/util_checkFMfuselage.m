%function fragment (see ForceAndMomentAero) debug and tuning
close all

alpvec = pi * (-1:0.001:1); %fuselage angle of attack

Fl_arr = NaN(2, length(alpvec));
Fn_arr = NaN(2, length(alpvec));

rho = 1.22;
Sfuseff = 0.4;
Cl_fus = 0.5;
Cd_fus = 0.01;
Cdmax = 1.1;
alpmax = 15*pi/180;

for n = 1:length(alpvec)
    Vb_uav = [cos(-alpvec(n)); sin(-alpvec(n))];
    
    Vorthog = sqrt(Vb_uav(2)^2);
    Ffus2D = FaeroMaero2D(-[Vb_uav(1); Vorthog], rho, Sfuseff, 0, 0, -Cl_fus, Cl_fus, Cd_fus, Cdmax, alpmax);
    Fl_arr(:,n) = [Ffus2D(1); Ffus2D(2) * Vb_uav(2)/(Vorthog + 1e-6)];
    
    Cnl = OrientMatrix2D(alpvec(n));
    Fn_arr(:,n) = Cnl * Fl_arr(:,n);
end



figure
subplot(211)
plot(180/pi*alpvec, Fn_arr(1,:), 'b', 180/pi*alpvec, Fn_arr(2,:), 'r'); grid on
subplot(212)
plot(180/pi*alpvec, 180/pi * (atan2(Fn_arr(2,:), Fn_arr(1,:)) - alpvec)); grid on

%check the direction of aero. force
hfig = figure;
for n = 1:length(alpvec)
    figure(hfig)
    plot([0 cos(alpvec(n))], [0 sin(alpvec(n))], 'g'); grid on; hold on
    plot([0 Fn_arr(1,n)], [0 Fn_arr(2,n)]); hold off
    axis(2*[-1 1 -1 1])
    drawnow
end










