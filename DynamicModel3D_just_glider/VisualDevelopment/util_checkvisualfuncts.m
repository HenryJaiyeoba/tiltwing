%to test visualization
close all

hax = axes;

Lcube = 10;

psi = 0;
theta = 0;
gamma = 0;
joy = vrjoystick(1);

visarrstruct = GenVisualArrays();

Ruav = [0; 5; 0];
dr = 0.4;

for n = 1:Inf
    joy_yaw = axis(joy,4);
    joy_pitch = axis(joy,2);
    joy_roll = axis(joy,1);
    phi_tilt = 2*pi/3 * 0.5*(1-axis(joy,3));
    
    psi = psi + 0.04*joy_yaw;
    theta = theta + 0.04*joy_pitch;
    gamma = gamma + 0.04*joy_roll;
    
    Ruav = Ruav + dr *[cos(theta)*cos(psi); sin(theta); cos(theta)*sin(psi)];
    
    ail_defl = joy_roll/2;
    elev_defl = joy_pitch/2;
    
    tic
    VisualizeUAV3D(Ruav, psi, theta, gamma, phi_tilt, ail_defl, elev_defl, visarrstruct, hax, Lcube);
    title(['aildefl = ' num2str(ail_defl*180/pi) ' deg.; elevdefl = ' num2str(elev_defl*180/pi) ' deg.'])
    drawnow
    toc
end