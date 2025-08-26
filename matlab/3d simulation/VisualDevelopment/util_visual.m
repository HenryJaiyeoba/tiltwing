%to prepare all visual details
close all

ang = pi * nthroot(((-10:10)/10), 3);
rad = exp((1*ang).^2);

x0vec = rad .* cos(ang);
y0vec = rad .* sin(ang);
B0max = max(abs(x0vec));
xrow = x0vec / B0max;
yrow = y0vec / B0max;

Broot = 0.5;
Btip = Broot/2;
Lhalfwing = 1.5;
Lailroot = 1; 
Lailtip = 1.3;
Bail = 0.7*Btip;

Rfuselage = 0.2;

%fuselage parameters
Lfus = 2.5;
Lxfusfront = 0.8;
fus_points = 10;
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
Xwingarr0right = [Broot*xrow; (Broot + (Btip-Broot)*Lailroot/Lhalfwing)*xrow; Bail*xrow; Bail*xrow; (Broot + (Btip-Broot)*Lailtip/Lhalfwing)*xrow; Btip*xrow] + Broot/5;
Ywingarr0right = [Broot*yrow; (Broot + (Btip-Broot)*Lailroot/Lhalfwing)*yrow; Bail*yrow; Bail*yrow; (Broot + (Btip-Broot)*Lailtip/Lhalfwing)*yrow; Btip*yrow];
Zwingarr0right = [0; Lailroot; Lailroot; Lailtip; Lailtip; Lhalfwing] * ones(size(xrow)) + Rfuselage;
Xwingarr0left = Xwingarr0right;
Ywingarr0left = -Ywingarr0right;
Zwingarr0left = -Zwingarr0right;

%AIRFOIL
%figure
%plot(xrow, yrow, '.-'); grid on; axis(5*[-1 1 -1 1])

%3D IMAGE
% figure
% surfl(Zwingarr0right, Xwingarr0right, Ywingarr0right); grid on; hold on;
% surfl(Zwingarr0left, Xwingarr0left, Ywingarr0left);
% surfl(Zfus0, Xfus0, Yfus0)
% axis(2*[-1 1 -1 1 -1 1])
% %shading interp
% colormap(copper)

%EXPERIMENT: draw in dynamics to see how fast it can be
gammavec = 0:0.01:3;
hfig = figure;
t0 = cputime;
for n = 1:length(gammavec)
    gamma = gammavec(n);
    Cnb = OrientMatrix(0.3, 0.2, gamma);
    
    Xfus = Cnb(1,1)*Xfus0 + Cnb(1,2)*Yfus0 + Cnb(1,3)*Zfus0;
    Yfus = Cnb(2,1)*Xfus0 + Cnb(2,2)*Yfus0 + Cnb(2,3)*Zfus0;
    Zfus = Cnb(3,1)*Xfus0 + Cnb(3,2)*Yfus0 + Cnb(3,3)*Zfus0;
    
    Xwingarrleft = Cnb(1,1)*Xwingarr0left + Cnb(1,2)*Ywingarr0left + Cnb(1,3)*Zwingarr0left;
    Ywingarrleft = Cnb(2,1)*Xwingarr0left + Cnb(2,2)*Ywingarr0left + Cnb(2,3)*Zwingarr0left;
    Zwingarrleft = Cnb(3,1)*Xwingarr0left + Cnb(3,2)*Ywingarr0left + Cnb(3,3)*Zwingarr0left;
    Xwingarrright = Cnb(1,1)*Xwingarr0right + Cnb(1,2)*Ywingarr0right + Cnb(1,3)*Zwingarr0right;
    Ywingarrright = Cnb(2,1)*Xwingarr0right + Cnb(2,2)*Ywingarr0right + Cnb(2,3)*Zwingarr0right;
    Zwingarrright = Cnb(3,1)*Xwingarr0right + Cnb(3,2)*Ywingarr0right + Cnb(3,3)*Zwingarr0right;
    
    figure(hfig)
    surf(Zfus, Xfus, Yfus); grid on; hold on;
    surf(Zwingarrleft, Xwingarrleft, Ywingarrleft);
    surf(Zwingarrright, Xwingarrright, Ywingarrright);
    title(['n = ' num2str(n) '; time = ' num2str(cputime-t0)])
    axis(2*[-1 1 -1 1 -1 1])
    %shading interp
    colormap(copper)
    drawnow
    hold off
    
end









