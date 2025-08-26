function A = CrossProdMatr(avec)
%matrix representing cross-product by vector avec
A=[0 -avec(3) avec(2); avec(3) 0 -avec(1); -avec(2) avec(1) 0];