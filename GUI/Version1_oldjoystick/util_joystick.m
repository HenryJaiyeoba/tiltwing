%just to see what is what in joystick
close all

joy = vrjoystick(1);

while 1
    controls = axis(joy, 1:5);
    plot(controls, 'o-'); grid on
    axis([-Inf Inf -1 1])
    drawnow
end