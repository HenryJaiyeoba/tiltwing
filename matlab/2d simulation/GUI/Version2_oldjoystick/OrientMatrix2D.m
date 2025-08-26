function Cnb = OrientMatrix2D(theta)
%orientation (rotation) matrix

c11 = cos(theta);
c12 = -sin(theta);
c21 = sin(theta);
c22 = cos(theta);

Cnb = [c11, c12; c21, c22];

