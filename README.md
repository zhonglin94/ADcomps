# ADcomps
ADcomps is a multi-phase compositional simulator which is able to model and evaluate the development performance of CO2 injection considering the adverse effect CO2 solubility. Finding a creative algorithm that correlates the phase behavior of oil-water(polar-nonpolar) system and solving non-linear flow equations through Newton iteration efficiently are the two most challenging parts. This simulator uses more rigorous fluid phase equilibrium theories(EOS+Wong-Sandler mixing rule+NRTL model) to consider the solubility of CO2 in the aqueous phase compared with leading commercial simulators. ADcomps followed the coding convention of MRST which is a valuable open source simulator embedding in MATLAB. This simulator is developed base on the frame of MRST 2017a and it should be run with the corresponding version.
# Installation
1. MRST 2017a(https://www.sintef.no/Projectweb/MRST/);
2. ADcomps (add the path of ADcomps into the modules of MRST);
# Try the ADcomps example
1. Run the startup in MRST;
2. Run the compositionalExample in ADcomps.


![Image text](https://github.com/zhonglin94/ADcomps/blob/master/image%20folder/SPE_1%20example.gif)
Fig. 1 CO2 solubility in water (mol%) during CO2 flooding calculated by ADcomps

![Image text](https://github.com/zhonglin94/ADcomps/blob/master/image%20folder/Simulation%20results%20of%20ADcomps.png)

Fig. 2 The performance of ADcomps compared with ECLIPSE where water is treated as an inert phase
# The framework of ADcomps
![Image text](https://github.com/zhonglin94/ADcomps/blob/master/image%20folder/ADcomps%20framework.png)
Note: ADcomps moudules are marked with red flags
# Notice
1. I upload the ECE as a zip file. It should be uncompressed before running the example;
2. It sometimes goes wrong because of the path. If the scripts fail to add the files automatically, please add them into the MATLAB environment manually. 
