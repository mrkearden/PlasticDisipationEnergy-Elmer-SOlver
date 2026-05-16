# PlasticDisipationEnergy-Elmer-SOlver  
PDESOlver.F90 not used anymore, all calculations are in UMAT  

Soler for Elmer for plastic energy disipation and heat increase  
Uses UMAT for perfectly plastic material model
Uses new solver to calculate plastic work, heat increases, and cummlative T increase  

Cuurent output is a csv file openable in a spreadsheet that hsa the results for each
element, npt, and ntens.  Example problem herein is a solid block with uniaxial  
stress in the Z direction,  so most other direction has small results, but Ntens 3 has 
verified results.  

