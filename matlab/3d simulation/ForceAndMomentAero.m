function [Fbvec, Mbvec, Fb_rail, Fb_lail, Fb_elev, Mb_rail, Mb_lail, Mb_elev, Mb_rwing, Mb_lwing] = ForceAndMomentAero(Vb_uav, omegab_uav, phiwing, phiaill, phiailr, phielev, rbvec_masscenter, AeroPars_struct)
%this function computes the aerodynamic force and momentum (with respect
%to mass center) in projection on RF "b"

%Vb_uav - UAV airspeed
%AeroPars_struct - all data describing UAV aerodynamic features should be stored here

rho = AeroPars_struct.rho; %air density
Shalfwing = AeroPars_struct.Shalfwing; %halfwing area 
MAC = AeroPars_struct.MAC; %wing MAC
rb_rhalfwingAC = AeroPars_struct.rb_rhalfwingAC; %location of the aerod. center of right halfwing (at "effective" distance from longitudal axis)
Sail = AeroPars_struct.Sail; %aileron area
rb_railAC = AeroPars_struct.rb_railAC; %location of the right aileron (aerod. center - assume does not change with phi_ail)
Selev = AeroPars_struct.Selev; %aileron area
rb_elevAC = AeroPars_struct.rb_elevAC; %location of the right aileron (aerod. center - assume does not change with phi_elev)
%airfoil params (cosider them the same for all aerod. surfaces)
alpmax = AeroPars_struct.alpmax; %stall angle
Cdmin = AeroPars_struct.Cdmin; %min. drag coef.
Cdmax = AeroPars_struct.Cdmax; %max. drag coef.
Clmax = AeroPars_struct.Clmax; %max. lift coefficient
Clmin = AeroPars_struct.Clmin; %min. lift coefficient (should be negative)
Cm = AeroPars_struct.Cm; %"nose down" moment coef.
%fuselage aerod. params (for very crude consideration)
Sfuseff = AeroPars_struct.Sfuseff; %"effective" fuselage area - to treat it as a wing to a certain extent (take axial symmetry into account)
rb_fusAC = AeroPars_struct.rb_fusAC; %REFINE!!!
Cl_fus = AeroPars_struct.Cl_fus; %REFINE!!!
Cd_fus = AeroPars_struct.Cd_fus; %REFINE!!!

%find all local airspeeds
rb_lhalfwingAC = [rb_rhalfwingAC(1); rb_rhalfwingAC(2); -rb_rhalfwingAC(3)]; 
rb_lailAC = [rb_railAC(1); rb_railAC(2); -rb_railAC(3)]; 

Vb_rwing = Vb_uav + cross(omegab_uav, rb_rhalfwingAC - rbvec_masscenter);
Vb_lwing = Vb_uav + cross(omegab_uav, rb_lhalfwingAC - rbvec_masscenter);
Vb_rail = Vb_uav + cross(omegab_uav, rb_railAC - rbvec_masscenter);
Vb_lail = Vb_uav + cross(omegab_uav, rb_lailAC - rbvec_masscenter);
Vb_elev = Vb_uav + cross(omegab_uav, rb_elevAC - rbvec_masscenter);
%airspeed in vertical cross-section (XY)
Vb_rwing2D = Vb_rwing(1:2); 
Vb_lwing2D = Vb_lwing(1:2); 
Vb_rail2D = Vb_rail(1:2); 
Vb_lail2D = Vb_lail(1:2); 
Vb_elev2D = Vb_elev(1:2);

Cbw = OrientMatrix2D(phiwing); %for the wing
Cba_l = OrientMatrix2D(phiwing+phiaill); %for the left aileron
Cba_r = OrientMatrix2D(phiwing+phiailr); %for the right aileron
Cbe = OrientMatrix2D(phielev); %for the elevator

%2D airspeeds of aerod. surfaces in their own local RFs
Vw_rwing2D = Cbw' * Vb_rwing2D;
Vw_lwing2D = Cbw' * Vb_lwing2D;
Var_rail2D = Cba_r' * Vb_rail2D;
Val_lail2D = Cba_l' * Vb_lail2D;
Ve_elev2D = Cbe' * Vb_elev2D;

%calculate 2D forces and moments
[Fw_rwing2D, Mrwing] = FaeroMaero2D(-Vw_rwing2D, rho, Shalfwing, MAC, Cm, Clmin, Clmax, Cdmin, Cdmax, alpmax);
[Fw_lwing2D, Mlwing] = FaeroMaero2D(-Vw_lwing2D, rho, Shalfwing, MAC, Cm, Clmin, Clmax, Cdmin, Cdmax, alpmax);
%ignore aerod. torque from steering surfaces
Far_rail2D = FaeroMaero2D(-Var_rail2D, rho, Sail, 0, 0, Clmin, Clmax, Cdmin, Cdmax, alpmax);
Fal_lail2D = FaeroMaero2D(-Val_lail2D, rho, Sail, 0, 0, Clmin, Clmax, Cdmin, Cdmax, alpmax);
Fe_elev2D = FaeroMaero2D(-Ve_elev2D, rho, Selev, 0, 0, Clmin, Clmax, Cdmin, Cdmax, alpmax);

%fuselage effect
Vorthog = sqrt(Vb_uav(2)^2 + Vb_uav(3)^2);
Ffus2D = FaeroMaero2D(-[Vb_uav(1); Vorthog], rho, Sfuseff, 0, 0, -Cl_fus, Cl_fus, Cd_fus, Cdmax, alpmax); %REFINE use of Cdmax and alpmax (however this needs CFD)
Fb_fus = [Ffus2D(1); Ffus2D(2) * Vb_uav(2)/(Vorthog + 1e-6); Ffus2D(2) * Vb_uav(3)/(Vorthog + 1e-6)];

%calculate all forces
Fb_rwing = [Cbw * Fw_rwing2D; 0];
Fb_lwing = [Cbw * Fw_lwing2D; 0];
Fb_rail = [Cba_r * Far_rail2D; 0];
Fb_lail = [Cba_l * Fal_lail2D; 0];
Fb_elev = [Cbe * Fe_elev2D; 0];

Fbvec = Fb_rwing + Fb_lwing + Fb_rail + Fb_lail + Fb_elev + Fb_fus;

%calculate all moments
Mb_rwingAC = [0; 0; Mrwing];
Mb_lwingAC = [0; 0; Mlwing];

Mb_rwing = cross(rb_rhalfwingAC - rbvec_masscenter, Fb_rwing);
Mb_lwing = cross(rb_lhalfwingAC - rbvec_masscenter, Fb_lwing);
Mb_rail = cross(rb_railAC - rbvec_masscenter, Fb_rail);
Mb_lail = cross(rb_lailAC - rbvec_masscenter, Fb_lail);
Mb_elev = cross(rb_elevAC - rbvec_masscenter, Fb_elev);
Mb_fus = cross(rb_fusAC - rbvec_masscenter, Fb_fus);

Mbvec = Mb_rwing + Mb_lwing + Mb_rail + Mb_lail + Mb_elev + Mb_fus + Mb_rwingAC + Mb_lwingAC;

































