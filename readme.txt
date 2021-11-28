
MM-FOD.R
The MM-FOD.R file contains the proposed MM-FOD algorithm for detecting functional outliers. 
In the bottom of the file, an illustrative example is shown to help the readers to better use the MM-FOD algorithm. Specifically, 
we show the input, the step by step calculation of the MM-FOD method, and the output of the algorithm.


simulation code for MM-FOD.rar
This rar file is for the simulation study in the paper "Data adaptive functional outlier detection: Analysis of the Paris bikesharing system data".
We compared the proposed MM-FOD method to 5 state-of-art methods, they are 
(1) the step-wise functional outlier detection (SFOD) method proposed by Yu et al. (2012) 
(2) the refined least trimmed functional scores (RLTFS) method proposed by Ren et al. (2017); 
(3) Graphical methods: the functional boxplot method with two most popular definitions of depth: the modified band depth (FB-MBD) proposed by L´opez-Pintado and
Romo (2009) and the total variation depth (FB-TVD) by Huang and Sun (2019); 
(4) the Archetypoid analysis (ADA) method proposed by Vinue and Epifanio (2021) ;
(5) the outliergram (OUG)approach by Arribas-Gil and Romo (2014).

Specifically, the file code_main.R contains the simulation setup, 
the file fun_simulation.R contains the implementation of the proposed MM-FOD algorithm and the five benchmark functional outlier detection methods,
the file allfuns.R contains the subfunctions needed for the above methods, which are directly called in the  fun_simulation.R file.
Futhermore, the "code for ada and oug method" includes the code for the ADA and OUG methods.