# Power System Security Enhancement

All the codes in this repository are related to the enhancement of power system security for real-time power system operations. The two main programs in the repository are as follows:
(1) Main_first_component.m: This program secures the power system aganist post-contingency cut-set saturation as well as critical branch overloads, during successive transmission outage scenarios. 
(2) Main_second_component.m: This program secures the power system only against post-contingency cut-set saturation, during successive transmission outage scenarios.
These main programs utilize several user defined functions/subroutines contained in the repository. MATLAB 2017 and Gurobi (version 9) were used for developing the algorithms. 

MATLAB Version:
The project can be run using version 2017a or higher.

Additional software related information:
1. MATPOWER 
2. GUROBI Optimizer
3. For better computational speed please comment out all the print statements in rundcpf.m and printpf.m

Theoretical foundations of these algorithms can be found in the following documents (available in the repository named "Papers").
1. R. Sen Biswas, A. Pal, T. Werho, and V. Vittal, "A Graph Theoretic Approach to Power System Vulnerability Identification," IEEE Transactions on Power Systems, vol. 36, no. 2, pp. 923-935, March 2021.
2. R. Sen Biswas, A. Pal, T. Werho, and V. Vittal, "Mitigation of Saturated Cut-sets During Multiple Outages to Enhance Power System Security," IEEE Transactions on Power Systems, April 2021.
3. R. Sen Biswas, A. Pal, T. Werho, and V. Vittal, "Fast Identification of Saturated Cut-sets using Graph Search Techniques," 2020 IEEE Power & Energy Society General Meeting (PESGM), pp. 1-5, Aug. 2020. 
4. R. Sen Biswas, "Power System Security Enhancement for Real-Time Operations During Multiple Outages using Network Science", PhD Dissertation, Arizona State University, May. 2021.

If you use the algorithms contained in this repository for your research and development, please cite the above mentioned papers. In case you have questions, please contact Reetam Sen Biswas at rsenbisw@asu.edu.


  
   
