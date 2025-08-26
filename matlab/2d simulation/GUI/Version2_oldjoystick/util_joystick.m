%just to see what is what in joystick
close all

joy = vrjoystick(1);

hfig = figure;

while 1
    controls = axis(joy);
    figure(hfig)
    plot(controls, 'o-'); grid on
    axis([-Inf Inf -1 1])
    drawnow
end