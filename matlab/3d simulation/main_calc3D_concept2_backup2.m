%computation of flight dynamics
%all values in SI

%the inertial reference frame (RF "i") OXYZ_i is similar to the navigational 
%RF PNUE: second AXIS OY is vertical

%NOTE: here we cotrol tiltwing but no thrust and no suck speed
%consideration

close all

%INCORPORATE JOYSTICK
joy = vrjoystick(1);
phiailmax = 25 * pi/180; %for ailerons
phielevmax = 15 * pi/180; %for elevator

%the majority of "physics" and computational parameters are here
for dummy = 1
%general consts
g_grav = 9.81;
Vb_air = [0;0;0];
%temporal parameters
dt = 1e-3; %time step in computation
Tcalc = 1000;

%on-the-ground params
theta_ground = 10 * pi/180;
Hground = 0.5;
    
%UAV cofig. parameters
%NOTE: we assume tiltwing axis to be the origin of RF "b" when setting
%positions of various keypoints

m_UAV = 45; %UAV mass
rbvec_masscenter = [-0.03; -0.15; 0]; %mass center
Lfus = 2; %fuselage length
Dfus = 0.45; %fuselage width
Lwing = 3.3; %total wingspan (for crude evaluation of inertia tensor)
maxchord = 0.5; %m, root chord (for crude evaluation of inertia tensor)
maxthickness = 0.06; %m, wing max. thickness (for crude evaluation of inertia tensor)
Jb_UAV = CalcFrameInTensor(m_UAV, Lfus, Dfus, Lwing, maxchord, maxthickness); %inertia tensor
Rprop = 0.5; %propeller radius
m_prop = 0.4; %propeller mass (to simulate limitations of change of rot. speed)

%aerodynamic parameters:
AeroPars_struct = struct;
AeroPars_struct.rho = 1.225; %air density
AeroPars_struct.Shalfwing = 0.56; %halfwing area 
AeroPars_struct.MAC = 0.4; %wing MAC
AeroPars_struct.rb_rhalfwingAC = [0; 0; 1.1]; %location of the aerod. center of right halfwing (at "effective" distance from longitudal axis)
AeroPars_struct.Sail = 0.04; %aileron area
AeroPars_struct.rb_railAC = [-0.15; 0; 1.5]; %location of the right aileron (aerod. center - assume does not change with phi_ail)
AeroPars_struct.Selev = 0.12; %elevator area
AeroPars_struct.rb_elevAC = [-1.5; -0.15; 0]; %location of the right aileron (aerod. center - assume does not change with phi_elev)
%airfoil params (cosider them the same for all aerod. surfaces)
AeroPars_struct.alpmax = 15*pi/180; %stall angle
AeroPars_struct.Cdmin = 0.02; %min. drag coef.
AeroPars_struct.Cdmax = 1.1; %max. drag coef.
AeroPars_struct.Clmax = 1.4; %max. lift coefficient
AeroPars_struct.Clmin = -1.2; %min. lift coefficient (should be negative)
AeroPars_struct.Cm = -0.03; %"nose down" moment coef.
%fuselage aerod. params (for very crude consideration)
AeroPars_struct.Sfuseff = 0.4; %"effective" fuselage area - to treat it as a wing to a certain extent (take axial symmetry into account)
AeroPars_struct.rb_fusAC = [-0.2; -0.15; 0]; %REFINE!!!
AeroPars_struct.Cl_fus = 0.5; %REFINE!!!
AeroPars_struct.Cd_fus = 0.01; %REFINE!!!

%[theta_set, Treq] = FindOptAlphaReqThrust(Vset, m_UAV, AeroPars_struct); %determine required attack angle and thrust

end

%INITIALIZATION:
%initial UAV parameters
psi0 = 0.2;
theta0 = 0; 
gamma0 = 0;
R0vec = [0; 100; 0];
V0vec = [20; 0; 0];
omega0b = [0; 0; 0]; %initial angular velocity (with respect to RF "n") in body RF coords

%visualization params
deltat_draw = 0.05; %time step in visualization
Lcube = 15;
Ntrajdraw = 20;
deltat_traj = 0.2;
Mdraw = round(deltat_draw/dt);
Mtraj = round(deltat_traj/dt);
visarrstruct = GenVisualArrays();
Rtrajarr = repmat(R0vec, 1, Ntrajdraw);
%prepare arrays and parameters for computation
tvec = 0:dt:Tcalc;
Ncalc = length(tvec);

%MAIN LOOP:
%arrays for the main loop
Rnarr = NaN(3,Ncalc);
Vnarr = NaN(3,Ncalc);
psiarr = NaN(1,Ncalc);
thetaarr = NaN(1,Ncalc);
gammaarr = NaN(1,Ncalc);
omegabnarr = NaN(3,Ncalc);
phiwingarr = NaN(1,Ncalc);
phiaillarr = NaN(1,Ncalc);
phiailrarr = NaN(1,Ncalc);
phielevarr = NaN(1,Ncalc);
%initialize
Rnarr(:,1) = R0vec;
Vnarr(:,1) = V0vec;
psiarr(1) = psi0;
thetaarr(1) = theta0;
gammaarr(1) = gamma0;
omegabnarr(:,1) = omega0b;
    %omega_lproparr(1) = omega_lprop0;
    %omega_rproparr(1) = omega_rprop0;
    
%to start the loop
Rn_1 = Rnarr(:,1);
Vn_1 = Vnarr(:,1);
Cnb_1 = OrientMatrix(psiarr(1), thetaarr(1), gammaarr(1));
omegab_1 = omegabnarr(:,1);

figure
hax = axes;

for n = 2:Ncalc
    %tic
    
    %take readings from joystick
    joy_pitch = axis(joy,2);
    joy_roll = axis(joy,1);
    phiwing = 2*pi/3 * 0.5*(1-axis(joy,3));
    phiailr = -joy_roll*phiailmax;
    phiaill = -phiailr;
    phielev = -joy_pitch*phielevmax;
    
       
    %PHYSICS
    Vb_1 = Cnb_1' * Vn_1;
    [Fb_aero, Mb_aero, Fb_rail, Fb_lail, Fb_elev, Mb_rail, Mb_lail, Mb_elev, Mb_rwing, Mb_lwing] = ForceAndMomentAero(Vb_1, omegab_1, phiwing, phiaill, phiailr, phielev, rbvec_masscenter, AeroPars_struct);
    
    
    Mb_net = Mb_aero;% + Mb_prop + Mb_gyro;
    Lb = Jb_UAV * omegab_1;
    omegab = omegab_1 + Jb_UAV \ (Mb_net - CrossProdMatr(omegab_1) * Lb) * dt;
    
    Cb_1b = CalcRotMatrix(omegab*dt);
        
    Cnb = Cnb_1 * Cb_1b;
    %orthogonalization
    Esym = 0.5*(Cnb*Cnb' - eye(3));
    Cnb = (eye(3) - Esym) * Cnb;
    
    Fn_aero = Cnb*Fb_aero;
    Fn_net = m_UAV*[0; -g_grav; 0] + Fn_aero;% + Fn_prop;
    Vn = Vn_1 + Fn_net/m_UAV * dt;
    Rn = Rn_1 + Vn*dt;
    
    %calculate orient. angles
    gamma = atan2(-Cnb(2,3), Cnb(2,2));
    theta = atan(Cnb(2,1)/sqrt(Cnb(2,2)^2+Cnb(2,3)^2));
    psi = atan2(Cnb(3,1), Cnb(1,1));
    
    %on the ground:
    if Rn(2) <= Hground
        Rn(2) = Hground;
        Vn(2) = max(Vn(2), 0);
        theta = theta_ground;
        gamma = 0;
        Cnb = OrientMatrix(psi, theta, gamma);
        omegab = zeros(3,1);
    end
    
    %TO FINALIZE ITERATION
    %write down to arrays
%     Rnarr(:,n) = Rn;
%     Vnarr(:,n) = Vn;
%     gammaarr(n) = atan2(-Cnb(2,3), Cnb(2,2));
%     thetaarr(n) = atan(Cnb(2,1)/sqrt(Cnb(2,2)^2+Cnb(2,3)^2));
%     psiarr(n) = atan2(Cnb(3,1), Cnb(1,1));
%     omegabnarr(:,n) = omegab;
%     
%     phi_lductarr(n) = phi_lduct;
%     phi_rductarr(n) = phi_rduct;
%     omega_lductarr(n) = omega_lduct;
%     omega_rductarr(n) = omega_rduct;
%     omega_lproparr(n) = omega_lprop;
%     omega_rproparr(n) = omega_rprop;
%     Pcons_larr(n) = Pcons_l;
%     Pcons_rarr(n) = Pcons_r;
    
    %update previous values
    Rn_1 = Rn;
    Vn_1 = Vn;
    omegab_1 = omegab;
    Cnb_1 = Cnb;
    %toc
    
    %VISUALIZE
    for dummy = 1
    %update trajectory points
    if ((n-1)/Mtraj - round((n-1)/Mtraj)) == 0
        Rtrajarr = [Rtrajarr(:, 2:end), Rn];
    end 
        
    if ((n-1)/Mdraw - round((n-1)/Mdraw)) == 0
        VisualizeUAV3D(Rn, psi, theta, gamma, phiwing, phiaill, -phielev, visarrstruct, hax, Lcube);
        %add speed visualization and some other features
        hold on
        
        plot3([Rtrajarr(3,:), Rn(3)], [Rtrajarr(1,:), Rn(1)], [Rtrajarr(2,:), Rn(2)], 'b.-'); %plot recent trajectory
        
        %plot "reference" points for better speed perception
%         Lcubestep = 0.7*Lcube;
%         xrefgr = Lcubestep*round(Rn(1)/Lcubestep); yrefgr = Lcubestep*round(Rn(2)/Lcubestep); zrefgr = Lcubestep*round(Rn(3)/Lcubestep); 
%         xrefgr = xrefgr + Lcubestep*(-1:1); yrefgr = yrefgr + Lcubestep*(-1:1); zrefgr = zrefgr + Lcubestep*(-1:1); 
%         [xrefarr, yrefarr, zrefarr] = ndgrid(xrefgr, yrefgr, zrefgr);
%         xrefarr = reshape(xrefarr, [], 1); yrefarr = reshape(yrefarr, [], 1); zrefarr = reshape(zrefarr, [], 1);
%         plot3(zrefarr, xrefarr, yrefarr, 'm*')
        
        plot3(Rn(3) + [0 Vn(3)], Rn(1) + [0 Vn(1)], Rn(2) + [0 Vn(2)], 'y', 'linewidth', 2); hold off
        title(['t = ' num2str(tvec(n)) ' sec; V = ' num2str(3.6*norm(Vn)) ' km/h'])
        drawnow
    end
    end
end

































