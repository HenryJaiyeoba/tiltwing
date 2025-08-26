function [] = VisualizeUAV3D(Ruav, psi, theta, gamma, phi_tilt, lail_defl, rail_defl, elev_defl, visarrstruct, hax, Lcube)
%dynamic visualization of UAV in 3D flight

%ail_defl, elev_defl - deflection angles of ailerons and elevator (positive when roll/pitch increase)

psicam = psi + pi;
xcamera = cos(psicam);
zcamera = sin(psicam);
ycamera = 0.5;

%to depict the ground
Hgroundgr = 20; %maximal altitude, at which we still show ground
altitude_coef = min(1, exp(-0.1*(Ruav(2)-0.75*Lcube)));
dXZground = Lcube/8;
xcenground = dXZground * round(Ruav(1) / dXZground);
zcenground = dXZground * round(Ruav(3) / dXZground);
xgroundticks = xcenground + altitude_coef*dXZground*(-4:4);
zgroundticks = zcenground + altitude_coef*dXZground*(-4:4);
Xground1 = [xgroundticks; xgroundticks]; Zground1 = Ruav(3) + repmat(altitude_coef*Lcube/2*[-1; 1], 1, length(xgroundticks)); Yground1 = max(zeros(size(Xground1)), Ruav(2)-0.75*Lcube);
Zground2 = [zgroundticks; zgroundticks]; Xground2 = Ruav(1) + repmat(altitude_coef*Lcube/2*[-1; 1], 1, length(zgroundticks)); Yground2 = max(zeros(size(Xground2)), Ruav(2)-0.75*Lcube);

%take params and arrays from structure
Xfus0 = visarrstruct.Xfus0;
Yfus0 = visarrstruct.Yfus0;
Zfus0 = visarrstruct.Zfus0;
Xwingarr00left = visarrstruct.Xwingarr00left;
Ywingarr00left = visarrstruct.Ywingarr00left;
Zwingarr00left = visarrstruct.Zwingarr00left;
Xwingarr00right = visarrstruct.Xwingarr00right;
Ywingarr00right = visarrstruct.Ywingarr00right;
Zwingarr00right = visarrstruct.Zwingarr00right;
Rfuselage = visarrstruct.Rfuselage;
Lxfusfront = visarrstruct.Lxfusfront;
Lfus = visarrstruct.Lfus;
Lhalfwing = visarrstruct.Lhalfwing;
Broot = visarrstruct.Broot;
Btip = visarrstruct.Btip;
Bail = visarrstruct.Bail;
Lailroot = visarrstruct.Lailroot;
Lailtip = visarrstruct.Lailtip;
Lhalftail = visarrstruct.Lhalftail;
Btail = visarrstruct.Btail;


Cnb = OrientMatrix(psi, theta, gamma);

Xfus = Cnb(1,1)*Xfus0 + Cnb(1,2)*Yfus0 + Cnb(1,3)*Zfus0 + Ruav(1);
Yfus = Cnb(2,1)*Xfus0 + Cnb(2,2)*Yfus0 + Cnb(2,3)*Zfus0 + Ruav(2);
Zfus = Cnb(3,1)*Xfus0 + Cnb(3,2)*Yfus0 + Cnb(3,3)*Zfus0 + Ruav(3);

%rotate wings according to tilt angle
Xwingarr0left = cos(phi_tilt)*Xwingarr00left - sin(phi_tilt)*Ywingarr00left;
Ywingarr0left = sin(phi_tilt)*Xwingarr00left + cos(phi_tilt)*Ywingarr00left + Rfuselage/2;
Zwingarr0left = Zwingarr00left;
Xwingarr0right = cos(phi_tilt)*Xwingarr00right - sin(phi_tilt)*Ywingarr00right;
Ywingarr0right = sin(phi_tilt)*Xwingarr00right + cos(phi_tilt)*Ywingarr00right + Rfuselage/2;
Zwingarr0right = Zwingarr00right;

Xwingarrleft = Cnb(1,1)*Xwingarr0left + Cnb(1,2)*Ywingarr0left + Cnb(1,3)*Zwingarr0left + Ruav(1);
Ywingarrleft = Cnb(2,1)*Xwingarr0left + Cnb(2,2)*Ywingarr0left + Cnb(2,3)*Zwingarr0left + Ruav(2);
Zwingarrleft = Cnb(3,1)*Xwingarr0left + Cnb(3,2)*Ywingarr0left + Cnb(3,3)*Zwingarr0left + Ruav(3);
Xwingarrright = Cnb(1,1)*Xwingarr0right + Cnb(1,2)*Ywingarr0right + Cnb(1,3)*Zwingarr0right + Ruav(1);
Ywingarrright = Cnb(2,1)*Xwingarr0right + Cnb(2,2)*Ywingarr0right + Cnb(2,3)*Zwingarr0right + Ruav(2);
Zwingarrright = Cnb(3,1)*Xwingarr0right + Cnb(3,2)*Ywingarr0right + Cnb(3,3)*Zwingarr0right + Ruav(3);

%generate ailerons
l1 = Broot + (Btip-Broot)*Lailroot/Lhalfwing - Bail;
l2 = Broot + (Btip-Broot)*Lailtip/Lhalfwing - Bail;
Xail00right = [0, -l1*cos(rail_defl); 0, -l2*cos(rail_defl)] + Broot/5 - Bail;
Yail00right = [0, -l1*sin(rail_defl); 0, -l2*sin(rail_defl)];
Zail0right = [Lailroot, Lailroot; Lailtip, Lailtip] + Rfuselage;
Xail00left = [0, -l1*cos(lail_defl); 0, -l2*cos(lail_defl)] + Broot/5 - Bail;
Yail00left = [0, -l1*sin(lail_defl); 0, -l2*sin(lail_defl)];
Zail0left = -Zail0right;
%rotate ailerons according to tilt angle
Xail0left = cos(phi_tilt)*Xail00left - sin(phi_tilt)*Yail00left;
Yail0left = sin(phi_tilt)*Xail00left + cos(phi_tilt)*Yail00left + Rfuselage/2;
Xail0right = cos(phi_tilt)*Xail00right - sin(phi_tilt)*Yail00right;
Yail0right = sin(phi_tilt)*Xail00right + cos(phi_tilt)*Yail00right + Rfuselage/2;

Xailleft = Cnb(1,1)*Xail0left + Cnb(1,2)*Yail0left + Cnb(1,3)*Zail0left + Ruav(1);
Yailleft = Cnb(2,1)*Xail0left + Cnb(2,2)*Yail0left + Cnb(2,3)*Zail0left + Ruav(2);
Zailleft = Cnb(3,1)*Xail0left + Cnb(3,2)*Yail0left + Cnb(3,3)*Zail0left + Ruav(3);
Xailright = Cnb(1,1)*Xail0right + Cnb(1,2)*Yail0right + Cnb(1,3)*Zail0right + Ruav(1);
Yailright = Cnb(2,1)*Xail0right + Cnb(2,2)*Yail0right + Cnb(2,3)*Zail0right + Ruav(2);
Zailright = Cnb(3,1)*Xail0right + Cnb(3,2)*Yail0right + Cnb(3,3)*Zail0right + Ruav(3);

%generate elevator
Xtail0 = [0, -Btail*cos(elev_defl); 0, -Btail*cos(elev_defl)] + Lxfusfront - Lfus + Btail/2;
Ytail0 = [0, -Btail*sin(elev_defl); 0, -Btail*sin(elev_defl)];
Ztail0 = [-Lhalftail, -Lhalftail; Lhalftail, Lhalftail];
Xtail = Cnb(1,1)*Xtail0 + Cnb(1,2)*Ytail0 + Cnb(1,3)*Ztail0 + Ruav(1);
Ytail = Cnb(2,1)*Xtail0 + Cnb(2,2)*Ytail0 + Cnb(2,3)*Ztail0 + Ruav(2);
Ztail = Cnb(3,1)*Xtail0 + Cnb(3,2)*Ytail0 + Cnb(3,3)*Ztail0 + Ruav(3);


%VISUALIZE
axes(hax)
surfl(Zfus, Xfus, Yfus); grid on; hold on;
surfl(Zwingarrleft, Xwingarrleft, Ywingarrleft);
surfl(Zwingarrright, Xwingarrright, Ywingarrright);
surf(Zailright, Xailright, Yailright);
surf(Zailleft, Xailleft, Yailleft);
surf(Ztail, Xtail, Ytail);
if Ruav(2)<Hgroundgr
    plot3(Zground1, Xground1, Yground1, 'g', Zground2, Xground2, Yground2, 'g')
    plot3(Ruav(3), Ruav(1), 0, 'ko', 'markersize', 12)
end

    %title(['n = ' num2str(n) '; time = ' num2str(cputime-t0)])
%reduce stretch/compression in visualization
psigr = mod(psi, pi/2); %from 0 to pi/2
hor_coef = 1 / (cos(psigr) + sin(psigr));
axis([Ruav(3) + Lcube/2*hor_coef*[-1, 1], Ruav(1) + Lcube/2*hor_coef*[-1, 1], Ruav(2) + Lcube/2*[-1.5, 0.5]])

%shading interp
colormap(copper)
set(gca, 'CameraPosition', [Ruav(3); Ruav(1); Ruav(2)] + 1e3*[zcamera; xcamera; ycamera])
drawnow
hold off