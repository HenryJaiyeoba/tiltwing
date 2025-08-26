function Cab = CalcRotMatrix(phi_b_ab)
%compute rotation matrix from finite turn vector

Phi_b_ab = CrossProdMatr(phi_b_ab); 
phi = norm(Phi_b_ab); %include protection against 0/0
if phi < 1e-9 %protection against 0/0 (not needed in practical cases)
    sinphi_phi = 1;
    one_cosphi_phi2 = 0.5;
else
    sinphi_phi = sin(phi)/phi;
    one_cosphi_phi2 = (1 - cos(phi))/phi^2;
end
Cab = eye(3) + sinphi_phi * Phi_b_ab + one_cosphi_phi2 * Phi_b_ab^2;

