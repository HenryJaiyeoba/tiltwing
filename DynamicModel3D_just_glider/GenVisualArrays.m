function visarrstruct = GenVisualArrays()
%just to cut down visualization time: one-time calculation of required
%arrays

visarrstruct = struct;

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
visarrstruct.Lhalftail = Lhalfwing/3;
visarrstruct.Btail = Btip;

%fuselage parameters
Rfuselage = 0.2;
Lfus = 2.5;
Lxfusfront = 0.8;
fus_points = 20;
p1 = 0.4;
p2 = 0.75;

%write requirted params into structure
visarrstruct.Rfuselage = Rfuselage;
visarrstruct.Lhalfwing = Lhalfwing;
visarrstruct.Broot = Broot;
visarrstruct.Btip = Btip;
visarrstruct.Bail = Bail;
visarrstruct.Lailroot = Lailroot;
visarrstruct.Lailtip = Lailtip;
visarrstruct.Lxfusfront = Lxfusfront;
visarrstruct.Lfus = Lfus;


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

visarrstruct.Xfus0 = Xfus0;
visarrstruct.Yfus0 = Yfus0;
visarrstruct.Zfus0 = Zfus0;

%generate wing
visarrstruct.Xwingarr00right = [Broot*xrow; (Broot + (Btip-Broot)*Lailroot/Lhalfwing)*xrow; Bail*xrow; Bail*xrow; (Broot + (Btip-Broot)*Lailtip/Lhalfwing)*xrow; Btip*xrow] + Broot/5;
visarrstruct.Ywingarr00right = [Broot*yrow; (Broot + (Btip-Broot)*Lailroot/Lhalfwing)*yrow; Bail*yrow; Bail*yrow; (Broot + (Btip-Broot)*Lailtip/Lhalfwing)*yrow; Btip*yrow];
visarrstruct.Zwingarr00right = [0; Lailroot; Lailroot; Lailtip; Lailtip; Lhalfwing] * ones(size(xrow)) + Rfuselage;
visarrstruct.Xwingarr00left = visarrstruct.Xwingarr00right;
visarrstruct.Ywingarr00left = visarrstruct.Ywingarr00right;
visarrstruct.Zwingarr00left = -visarrstruct.Zwingarr00right;



