%this is a 2D model to analyze the transition mode of UAV
%! here only pitch control in hover mode
close all

%temporal parameters
dt = 1e-3; %time step in computation
Tcalc = 50;

%UAV parameters
m_UAV = 40; %UAV mass
LcharUAV = 0.5; %characteristic UAV length (determines inertia moment) 
Swing = 1.2;
Sstab = Swing/10;
MACwing = 0.4;
MACstab = 0.12;
Dprop = 1; %propeller diameter
coef_downwash = 0.7; %the portion of the wing area downwashed by propeller stream
alp_props = -6 * pi/180; %set angle of the motors (and propellers) with respect to wing chord
rb2D_acw = [0.02; 0.2]; %wing aero. center with respect to mass center
rb2D_acs = [-1.5; 0]; %wing aero. center with respect to mass center
%aerodynamic (consider the wing and stabiliter have the same airfoil)
Clmin = -1.2;
Clmax = 1.4;
Cdmin = 0.02;
Cdmax = 1.1;
Cm = -0.03;
alpmax = 15*pi/180; %stall angle
%servos
omegaservomax_wing = 60 * pi/180;
Mservomax_wing = 35; %maximal (rated) servo torque, N*m


%initial parameters
R0 = [0; 0.5];
V0 = [1; 0];
theta0 = 0 * pi/180;
omega0 = -20* pi/180;

%general parameters
g_grav = 9.81;
rho = 1.22;

%visualization parameters
Xfront = 0.6;
deltat_draw = 0.02; %time step in visualization
Lhalf_trajgr = 10; %box halfsize for trajectory picture

%prepare arrays and parameters for computation
tvec = 0:dt:Tcalc;
Ncalc = length(tvec);
I_UAV = m_UAV * LcharUAV^2; %inertia moment of the whole aircraft
Itiltwing = m_UAV/4 * (MACwing/2)^2; %inertia moment of the tiltwing
Sprop = 0.25*pi*Dprop^2;
phiwing0 = atan2(rb2D_acw(2), rb2D_acw(1)) - alp_props;
thetaset = pi/2 - atan2(rb2D_acw(2), rb2D_acw(1));

%prepare arrays and parameters for visualization
Mdraw = round(deltat_draw/dt);
Rbfusdraw = [rb2D_acs(1), Xfront; 0, 0];
[xairfoil, yairfoil] = GenAirfoilGrid();
xairfoil = xairfoil - 0.25; %to center the airfoil in aero. center
xairfoil = [xairfoil, 0];
yairfoil = [yairfoil, 0];
points_af = length(xairfoil);
Rwingdraw = 2 * MACwing * [xairfoil; yairfoil]; %numberic coef. for better visualization
Rstabdraw = 2 * MACstab * [xairfoil; yairfoil]; %numberic coef. for better visualization

%INITIALIZE
Rn_1 = R0;
Vn_1 = V0;
theta_1 = theta0;
omega_1 = omega0;
domegadt_cur = 0;
phiwing = 90 * pi/180;
dphiwingdt = 0;
Rrefgr = Rn_1; %for graphics

%MAIN LOOP
hfig = figure;
for n = 2:Ncalc
    % UAV dynamics
    Tprop = 1.01 * m_UAV*g_grav; %this is temporary - must be controlled
        %phiwing = 0 * pi/180; %this is temporary - must be controlled 
    phistab = -25 * pi/180; %this is temporary - must be controlled
    
    %tiltwing hover control (with real-world limitations)
    phiwingset = TiltwingHoverContoller(theta_1 - thetaset,  omega_1, phiwing0, norm(rb2D_acw), LcharUAV);
    dphiwingdt_set = (phiwingset - phiwing)/dt;
    [dphiwingdt, ~] = ConfineOmegaServo(dphiwingdt_set, dphiwingdt, dt, domegadt_cur, omegaservomax_wing, Mservomax_wing, Itiltwing);
    phiwing = phiwing + dt * dphiwingdt;
    
    Cnb_2D = OrientMatrix2D(theta_1);
    Cbw_2D = OrientMatrix2D(phiwing); %from wing to body
    Cbs_2D = OrientMatrix2D(phistab); %from elevator to body
    
    %thrust
    Fb_prop = Tprop * [cos(phiwing+alp_props); sin(phiwing+alp_props)];
    Fn_prop = Cnb_2D * Fb_prop;
    Mprop = rb2D_acw(1)*Fb_prop(2) - rb2D_acw(2)*Fb_prop(1);
    
    %calculate the increment of airspeed around the wing due to props
    %influence
    Vsuck = sqrt(0.5*Tprop / (rho * Sprop)); %"suck" speed - at this speed propeller "throws" air
    Vn_suck = -Vsuck * Fn_prop/norm(Fn_prop);
    
    %aerodynamic
    rw2D_acw = Cbw_2D' * rb2D_acw;
    Vw_air_wing = (Cnb_2D*Cbw_2D)' * (-Vn_1 + coef_downwash*Vn_suck); %airspeed for the wing is assumed the same as for UAV CG
    rs2D_acs = Cbs_2D' * rb2D_acs;
    rn2D_acs = Cnb_2D * rb2D_acs;
    Vs_air_stab = (Cnb_2D*Cbs_2D)' * -(Vn_1 + omega_1 * [-rn2D_acs(2); rn2D_acs(1)]); %unlike the wing, for the tail we consider the effect of rotary motion
    [Fw_wing, Mwing, ~, ~] = FaeroMaero(Vw_air_wing, rw2D_acw, rho, Swing, MACwing, Cm, Clmin, Clmax, Cdmin, Cdmax, alpmax);
    [Fs_stab, Mstab, ~, ~] = FaeroMaero(Vs_air_stab, rs2D_acs, rho, Sstab, MACstab, Cm, Clmin, Clmax, Cdmin, Cdmax, alpmax);
    Fn_wing = Cnb_2D*Cbw_2D * Fw_wing;
    Fn_stab = Cnb_2D*Cbs_2D * Fs_stab;
    
    %integration
    Vn = Vn_1 + dt * ((Fn_wing + Fn_stab + Fn_prop)/m_UAV + [0; -g_grav]);
    domegadt_cur = (Mwing + Mstab + Mprop) / I_UAV;
    omega = omega_1 + dt * domegadt_cur;
    Rn = Rn_1 + dt*Vn;
    theta = theta_1 + dt*omega;
    
    %update
    Rn_1 = Rn;
    Vn_1 = Vn;
    theta_1 = theta;
    omega_1 = omega;
    
    %visualization
    if ((n-1)/Mdraw - round((n-1)/Mdraw)) == 0
        Rnfusdraw = Cnb_2D * Rbfusdraw + repmat(Rn, 1, 2);
        Rnwingdraw = Cnb_2D * (Cbw_2D * Rwingdraw + repmat(rb2D_acw, 1, points_af)) + repmat(Rn, 1, points_af);
        Rnstabdraw = Cnb_2D * (Cbs_2D * Rstabdraw + repmat(rb2D_acs, 1, points_af)) + repmat(Rn, 1, points_af);
        
        if n>(Mdraw+1)
            delete(h0); delete(h1); delete(h2); delete(h3); delete(h4)
        end
        figure(hfig)
        h0 = plot(Rnfusdraw(1,:), Rnfusdraw(2,:), 'k', 'linewidth', 2); grid on; hold on
        h1 = plot(Rnwingdraw(1,:), Rnwingdraw(2,:), 'r', 'linewidth', 1);
        h2 = plot(Rnstabdraw(1,:), Rnstabdraw(2,:), 'r', 'linewidth', 1);
            
            debcoef1 = 1 / (m_UAV*g_grav);
            debcoef2 = 1;
            rn2D_acw = Cnb_2D * rb2D_acw;
            h3 = plot(Rn(1)+rn2D_acw(1)+[0 debcoef1*Fn_wing(1)], rn2D_acw(2)+Rn(2)+[0 debcoef1*Fn_wing(2)], 'g.--', 'linewidth', 1);
            h4 = plot(Rn(1)+[0 debcoef2*Vn(1)], Rn(2)+[0 debcoef2*Vn(2)], 'b.--', 'linewidth', 1);
            
        axis(Lhalf_trajgr*[-1, 1, -1, 1] + [Rrefgr(1), Rrefgr(1), Rrefgr(2), Rrefgr(2)]);
        if max(abs(Rn - Rrefgr)) >= (0.8*Lhalf_trajgr)
           Rrefgr = Rn + 0.7*(Rn - Rrefgr);
        end
        title(['Time = ' num2str(tvec(n)) ' s; UAV speed = ' num2str(norm(Vn)) ' m/s'])
        drawnow
    end
    
end













