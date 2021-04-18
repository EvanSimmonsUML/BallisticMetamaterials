function epsSub=subEps(lam0)

%Permittivity of substrate used in manuscript and experiments. 
A_InP = 6.255;
B1_InP = 2.316;
C1_InP = 0.6263^2;
B2_InP = 2.765;
C2_InP = 32.935^2;

epsSub = 1 +  A_InP + (B1_InP.*lam0.^2)/(lam0.^2 - C1_InP)+...
    (B2_InP*lam0.^2)/(lam0.^2 - C2_InP);
end
