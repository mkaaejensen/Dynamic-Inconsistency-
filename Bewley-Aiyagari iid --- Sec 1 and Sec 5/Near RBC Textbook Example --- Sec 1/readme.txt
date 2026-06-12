General notes: 

- The tsa CMA-ES implementation is recommended in "near RBC" (or "near deterministic") models
- Originally considered for teaching purposes (the near RBC version of the Aiyagari model is commonly introduced to students and being able to present the quasi-hyperbolic version is obviously useful)
- aiyagari_Sec1_First_best.m computes first best/DC solution and compares with EGP (endogenous grid point) algorithm
- aiyagari_Sec1_PanelA.m computes the first panel in the Introduction
- aiyagari_Sec1_PanelB.m computes the second panel in the Introduction with the saving mandate derived from the console output from aiyagari_Sec1_PanelA.m

- As in optimal growth models, restart is often useful to speed up convergence (see relevant readme for details)
- If tsa fails to reach an acceptable loss, restart (cancel and then loadoutput=1) with default (change algorithm   = 'tsa' to algorithm   = 'default')