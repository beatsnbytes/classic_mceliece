This repository contains a hardware/software co-design acceleration of Classic McEliece Key Encapsulation Mechanism (CM KEM).
It is based on the 3rd Round submission of the CM KEM to the NIST Post-Quantum Standardization contest.
The baseline implementation for the different security levels can be found in the "mceliece-20201010/baseline/mceliece<>" branches. 
The master branch contains the baseline\mceliece348864 (Level 1) software implementation. 
The baseline code is the C software implementation of CM KEM with the addition of custom time measurement functions for different parts of the application.
Baseline software implementations have been augmented with custom time measurement functions for different parts of the application.

To compile and run the CM KEM application you should use the build file provided.
Specific flags for compilation are ARCH=<arm, x86> (default is x86), OPT=<number> where you define the gcc optimization level, VEC=<y, n> where the user defines if the coade auto-vectorization feature of gcc will be used and KAT=<number> where the user defines how many consecutive Known Answer Tests (KATs) will be executed.
An example compilation looks like this "make kat OPT=3 VEC=y KAT=10" where we specify compilation of the CM KEM Known Answer Test KAT application for -O3 gcc optimization level with code-autovectorization included and for 10 consecutive KAT tests. 
After we can execute the binray "./kat" and see the timing measurements at the standard output.

Our code is tested in Ubuntu 18.04. The user needs to install the libssl-dev package for the compilation to succeed.

 
