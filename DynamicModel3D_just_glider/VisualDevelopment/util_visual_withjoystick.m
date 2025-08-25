%to prepare all visual details; commands from joystick read, also we check
%how camera position shift works
close all

ang = pi * nthroot(((-10:10)/10), 3);
rad = exp((1*ang).^2);

x0vec = rad .* cos(ang);
y0vec = rad .* sin(ang);
B0max = max(abs(x0vec));
xrow = x0vec / B0max;
yrow = y0vec / B0max;

%wing parameters
Broot = 0.5;
Btip = Broot/2;
Lhalfwing = 1.5;
Lailroot = 0.8; 
Lailtip = 1.3;
Bail = 0.7*Btip;

%tail parameters
Lhalftail = Lhalfwing/3;
Btail = Btip;

%fuselage parameters
Rfuselage = 0.2;
Lfus = 2.5;
Lxfusfront = 0.8;
fus_points = 20;
p1 = 0.4;
p2 = 0.75;

%generate fuselage
pvec = (0:fus_points)/fus_points;
xfusvec = Lxfusfront - Lfus + Lfus*pvec;
r0fusvec = NaN(1, length(pvec));
for n = 1:length(pvec)
    pcur = pvec(n);
    if pcur > p2
        r0fusvec(n) = sqrt( 1 - ((pcur-p2)/(1-p2))^2 );
    elseif pcur < p1
        r0fusvec(n) = 1 - ((pcur-p1)/p1)^2;
    else 
        r0fusvec(n) = 1;
    end
end
phifus = pi*(-1:0.2:1);
Nphi = length(phifus);
Xfus0 = NaN(length(pvec), Nphi);
Yfus0 = NaN(length(pvec), Nphi);
Zfus0 = NaN(length(pvec), Nphi);
for n = 1:length(pvec)
    Xfus0(n,:) = xfusvec(n);
    Yfus0(n,:) = Rfuselage * r0fusvec(n) .* sin(phifus);
    Zfus0(n,:) = Rfuselage * r0fusvec(n) .* cos(phifus);
end

%generate wing
Xwingarr00right = [Broot*xrow; (Broot + (Btip-Broot)*Lailroot/Lhalfwing)*xrow; Bail*xrow; Bail*xrow; (Broot + (Btip-Broot)*Lailtip/Lhalfwing)*xrow; Btip*xrow] + Broot/5;
Ywingarr00right = [Broot*yrow; (Broot + (Btip-Broot)*Lailroot/Lhalfwing)*yrow; Bail*yrow; Bail*yrow; (Broot + (Btip-Broot)*Lailtip/Lhalfwing)*yrow; Btip*yrow];
Zwingarr00right = [0; Lailroot; Lailroot; Lailtip; Lailtip; Lhalfwing] * ones(size(xrow)) + Rfuselage;
Xwingarr00left = Xwingarr00right;
Ywingarr00left = Ywingarr00right;
Zwingarr00left = -Zwingarr00right;

%MAIN LOOP
hfig = figure;
t0 = cputime;
psi = 0;
theta = 0;
gamma = 0;

joy = vrjoystick(1);

for n = 1:Inf
    
    joy_yaw = axis(joy,4);
    joy_pitch = axis(joy,2);
    joy_roll = axis(joy,1);
    phi_tilt = 2*pi/3 * 0.5*(1-axis(joy,3));
    
    psi = psi + 0.04*joy_yaw;
    theta = theta + 0.04*joy_pitch;
    gamma = gamma + 0.04*joy_roll;
    
    psicam = psi + pi;
    xcamera = cos(psicam);
    zcamera = sin(psicam);
    ycamera = 0.5;
    
    Cnb = OrientMatrix(psi, theta, gamma);
    
    Xfus = Cnb(1,1)*Xfus0 + Cnb(1,2)*Yfus0 + Cnb(1,3)*Zfus0;
    Yfus = Cnb(2,1)*Xfus0 + Cnb(2,2)*Yfus0 + Cnb(2,3)*Zfus0;
    Zfus = Cnb(3,1)*Xfus0 + Cnb(3,2)*Yfus0 + Cnb(3,3)*Zfus0;
    
    %rotate wings according to tilt angle
    Xwingarr0left = cos(phi_tilt)*Xwingarr00left - sin(phi_tilt)*Ywingarr00left;
    Ywingarr0left = sin(phi_tilt)*Xwingarr00left + cos(phi_tilt)*Ywingarr00left + Rfuselage/2;
    Zwingarr0left = Zwingarr00left;
    Xwingarr0right = cos(phi_tilt)*Xwingarr00right - sin(phi_tilt)*Ywingarr00right;
    Ywingarr0right = sin(phi_tilt)*Xwingarr00right + cos(phi_tilt)*Ywingarr00right + Rfuselage/2;
    Zwingarr0right = Zwingarr00right;
    
    Xwingarrleft = Cnb(1,1)*Xwingarr0left + Cnb(1,2)*Ywingarr0left + Cnb(1,3)*Zwingarr0left;
    Ywingarrleft = Cnb(2,1)*Xwingarr0left + Cnb(2,2)*Ywingarr0left + Cnb(2,3)*Zwingarr0left;
    Zwingarrleft = Cnb(3,1)*Xwingarr0left + Cnb(3,2)*Ywingarr0left + Cnb(3,3)*Zwingarr0left;
    Xwingarrright = Cnb(1,1)*Xwingarr0right + Cnb(1,2)*Ywingarr0right + Cnb(1,3)*Zwingarr0right;
    Ywingarrright = Cnb(2,1)*Xwingarr0right + Cnb(2,2)*Ywingarr0right + Cnb(2,3)*Zwingarr0right;
    Zwingarrright = Cnb(3,1)*Xwingarr0right + Cnb(3,2)*Ywingarr0right + Cnb(3,3)*Zwingarr0right;
    
    %generate ailerons
    l1 = Broot + (Btip-Broot)*Lailroot/Lhalfwing - Bail;
    l2 = Broot + (Btip-Broot)*Lailtip/Lhalfwing - Bail;
    Xail00right = [0, -l1*cos(joy_roll/2); 0, -l2*cos(joy_roll/2)] + Broot/5 - Bail;
    Yail00right = [0, l1*sin(joy_roll/2); 0, l2*sin(joy_roll/2)];
    Zail0right = [Lailroot, Lailroot; Lailtip, Lailtip] + Rfuselage;
    Xail00left = Xail00right;
    Yail00left = -Yail00right;
    Zail0left = -Zail0right;
    %rotate ailerons according to tilt angle
    Xail0left = cos(phi_tilt)*Xail00left - sin(phi_tilt)*Yail00left;
    Yail0left = sin(phi_tilt)*Xail00left + cos(phi_tilt)*Yail00left + Rfuselage/2;
    Xail0right = cos(phi_tilt)*Xail00right - sin(phi_tilt)*Yail00right;
    Yail0right = sin(phi_tilt)*Xail00right + cos(phi_tilt)*Yail00right + Rfuselage/2;
    
    Xailleft = Cnb(1,1)*Xail0left + Cnb(1,2)*Yail0left + Cnb(1,3)*Zail0left;
    Yailleft = Cnb(2,1)*Xail0left + Cnb(2,2)*Yail0left + Cnb(2,3)*Zail0left;
    Zailleft = Cnb(3,1)*Xail0left + Cnb(3,2)*Yail0left + Cnb(3,3)*Zail0left;
    Xailright = Cnb(1,1)*Xail0right + Cnb(1,2)*Yail0right + Cnb(1,3)*Zail0right;
    Yailright = Cnb(2,1)*Xail0right + Cnb(2,2)*Yail0right + Cnb(2,3)*Zail0right;
    Zailright = Cnb(3,1)*Xail0right + Cnb(3,2)*Yail0right + Cnb(3,3)*Zail0right;
    
    %generate elevator
    Xtail0 = [0, -Btail*cos(joy_pitch/2); 0, -Btail*cos(joy_pitch/2)] + Lxfusfront - Lfus + Btail/2;
    Ytail0 = [0, Btail*sin(joy_pitch/2); 0, Btail*sin(joy_pitch/2)];
    Ztail0 = [-Lhalftail, -Lhalftail; Lhalftail, Lhalftail];
    Xtail = Cnb(1,1)*Xtail0 + Cnb(1,2)*Ytail0 + Cnb(1,3)*Ztail0;
    Ytail = Cnb(2,1)*Xtail0 + Cnb(2,2)*Ytail0 + Cnb(2,3)*Ztail0;
    Ztail = Cnb(3,1)*Xtail0 + Cnb(3,2)*Ytail0 + Cnb(3,3)*Ztail0;
    
    
    %VISUALIZE
    figure(hfig)
    surfl(Zfus, Xfus, Yfus); grid on; hold on;
    surfl(Zwingarrleft, Xwingarrleft, Ywingarrleft);
    surfl(Zwingarrright, Xwingarrright, Ywingarrright);
    surf(Zailright, Xailright, Yailright);
    surf(Zailleft, Xailleft, Yailleft);
    surf(Ztail, Xtail, Ytail);
        
    title(['n = ' num2str(n) '; time = ' num2str(cputime-t0)])
    %reduce stretch/compression in visualization
    psigr = mod(psi, pi/2); %from 0 to pi/2
    axis(4*[[-1 1 -1 1] / (cos(psigr) + sin(psigr)), -1, 1])
       
    %shading interp
    colormap(copper)
    set(gca, 'CameraPosition', 1e3*[zcamera xcamera ycamera])
    drawnow
    hold off
    
end









