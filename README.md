##WilsonBoot

The files exclude theories based on different sets of crossing equations.  
_isingDemo.py_ reproduces 3D results as a sanity check  
_test_Identical_1D.py_ checks only a single crossing equation and reproduces known 1D results  
_test\_O(N)\_1D.py_ uses 3 crossing equations but gets the same bounds  
_test\_mixedZ2\_1D.py_ uses 5 equations, effectively copying the 3D ising bootstrap. It gives improved bounds  
_wilsonboot.py_ is where I check individual points using the full system of 9 equations  
_scanfullsystem.py_ does a bisection to find and compare bounds in mixedZ2 and fullsystem
