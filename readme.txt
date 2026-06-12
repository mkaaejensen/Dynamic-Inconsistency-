This repository contains the codebase for Time-Consistent Saving and the Distribution of Wealth 

- Detailed documentation of all scripts can be found in documentation.pdf
- Section 8 and Section 9 in documentation.pdf map different models and the figures in the paper to the appropriate folders
- A suggested roadmap for becoming familiar with the code is
	* Begin with optimal_growth.m in the folder optimal growth --- deterministic --- Section 3
		This computes Panel B, Figure 4 in the paper and provides verification of the basic algorithm etc by comparing with the (only) known 			analytical solution
	* Also in optimal_growth.m, try changing delta or beta to see how Markov policy is affected and gain familiarity. Optionally, run the additional
		scripts in the folder which provide context and generalizations (e.g. the ...sigma2.m script computes a case with jumps)
	* For the State Space Recursion algorithm, run the scripts in the root of the folder determ cons sav --- Section 4. These models are computationally
		cheap and allow for comparison of a variety of different parameter configurations.
	* The "textbook near-RBC" example in the paper's introduction (folder Bewley-Aiyagari iid --- Sec 1 and Sec 5/Near RBC Textbook Example --- Sec 1)
		is a good place to start for stochastic models
	* The rest of the code can be tested side by side with reading the paper (the matlab file names are mentioned in figure footnotes and it should
		be easy to find the appropriate folder --- if not, consult Section 9 in documentation.m


All comments and corrections are always very welcome

Martin K Jensen, June 2026 