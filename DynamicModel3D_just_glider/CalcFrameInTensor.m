function Jb = CalcFrameInTensor(m_total, Lfus, Dfus, Lwing, maxchord, maxthickness)
%this is temporary funct.: roghly estimate UAV inertia tensor (diagonal)
%in RF "b"; "fus" - fuselage
%it is assumed density is uniform; ducts are ignored 

vol_fus = 4/3*pi * 1/8*Lfus*Dfus^2;
vol_wing = 4/3*pi * 1/8*Lwing*maxchord*maxthickness;

m_fus = m_total * vol_fus / (vol_fus + vol_wing);
m_wing = m_total * vol_fus / (vol_fus + vol_wing);

Ix_fus = 2/5*m_fus * 1/4*Dfus^2;
Iyz_fus = 1/5*m_fus * 1/4*(Dfus^2 + Lfus^2);
Jb_fus = diag([Ix_fus, Iyz_fus, Iyz_fus]);

Ix_wing = 1/5*m_wing * 1/4*(Lwing^2 + maxthickness^2);
Iy_wing = 1/5*m_wing * 1/4*(Lwing^2 + maxchord^2);
Iz_wing = 1/5*m_wing * 1/4*(maxchord^2 + maxthickness^2);
Jb_wing = diag([Ix_wing, Iy_wing, Iz_wing]);

Jb = Jb_fus + Jb_wing;

%visualize the imitated body:
theta = pi*(-1:0.1:1);
phi = pi/2*(-1:0.1:1);
[theta, phi] = meshgrid(theta, phi);

xfus = Lfus/2 * cos(theta) .* cos(phi);
yfus = Dfus/2 * sin(phi);
zfus = Dfus/2 * sin(theta) .* cos(phi);
xwing = maxchord/2 * cos(theta) .* cos(phi);
ywing = maxthickness/2 * sin(phi);
zwing = Lwing/2 * sin(theta) .* cos(phi);

% figure
% surfl(zfus, xfus, yfus); hold on
% surfl(zwing, xwing, ywing); hold on
% xlabel('Z'); ylabel('X'); zlabel('Y');
% axis(Lwing*[-1 1 -1 1 -1 1])
% shading interp; colormap(gray)
% title('UAV body model for inertia tensor computation')




