//parameters
pz(1) = 0.005    //k_fqdoR»     constant of formation of the qdoR protein
pz(2)= 0.012    //k_fgfp»      constant of formation of GFProtein(mol/min-1)
pz(3) = 0.1    //K_dR»     First order rate constant for the degradation of repressor (min-1) 
pz(4) = 0.00035    //K_dgfp»        First order rate constant for the degradation of GFProtein (min-1)
pz(5) = 100D-9    //[Q]»        concentration of free quercetin in the cell (mol L-1)
pz(6) = 12D-9    //[R_g]_TOT»      Total concentration of Repressor gene in the cell (mol L-1)
pz(7)= 12D-9    //[BM_g]_TOT»     Total concentration of Biomarker gene in the cell (mol L-1)
pz(8) = 0.7    //K_dQ»      Dissociation constant for the “Q~R” complex
pz(9) = 0.04    //K_dNAR»      Dissociation constant for the “Promoter~R” complex
pz(10)= 0.8    //K_dBMg»     Dissociation constant for the "Biomarker~R" complex
pz(11) = 0.001    //mi»       Constant for volume increase
[VrSTR]= [0.1; 0.15; 0.16; 0.24; 0.33; 0.45]*0.2    //vr = velocity of production of repressor (mol repressor L-1 min-1)
