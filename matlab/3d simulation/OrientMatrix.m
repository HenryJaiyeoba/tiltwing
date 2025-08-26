function Cnb = OrientMatrix(psi, theta, gamma)
%orientation (rotation) matrix

c11 = cos(theta) * cos(psi);
c12 = -cos(gamma) * sin(theta) * cos(psi) - sin(gamma) * sin(psi);
c13 = sin(gamma) * sin(theta) * cos(psi) - cos(gamma) * sin(psi);
c21 = sin(theta);
c22 = cos(gamma) * cos(theta);
c23 = -sin(gamma) * cos(theta);
c31 = cos(theta) * sin(psi);
c32 = -cos(gamma) * sin(theta) * sin(psi) + sin(gamma) * cos(psi);
c33 = sin(gamma) * sin(theta) * sin(psi) + cos(gamma) * cos(psi);

Cnb = [c11, c12, c13; c21, c22, c23; c31, c32, c33];

