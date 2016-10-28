This is mc++, a program mainly used to simulate XY systems in 2D with dipolar interaction. This is the basis for the PhD of Dominik Schildknecht <dominik.schildknecht@psi.ch>
The program is not part of the ALPS project (http://alps.comp-phys.org/), however is designed to be linked against it. However for anything surrounding the main simulation as input file generation and evaluation of the raw data alpspython is used. The software is therefore under the ALPS licencse v1.1, distributed with this program.

The papers needed to cite to use this program are displayed in the beginning of a run. Please note that it is mandatory to cite these if the results are published.

Installation Guide:
    - Download the ALPS library (http://alps.comp-phys.org/) and install it according to the tutorial on the ALPS page
    - Set the environment variable ALPS_DIR to where ever you installed the library to (e.g. /opt/alps/share/alps, for this you need to do the make install with root privileges)
        sudo vim /etc/environment
    - (If the structure factor is heavily used) Download the FFTW library and install it (on Ubuntu easier:)
        sudo apt-get install libfftw3-dev 
    - Clone the mc++ program
        git clone https://github.com/domischi/mcpp.git
    - Change into mc++ folder
        cd mcpp
    - Go into the build directory and use cmake & make
        cd build; cmake ..; make
    - Add the two links to the alps programs and to mc++ to the PATH variable
        sudo vim /etc/environment
        (PATH="...:/opt/alps/bin:<mc++-install-directory>/bin")

User Guide:
    Follow mainly the example. If so you are ensured to generate files that are compatible with this program. You can either start the program directly in python as given in the examples, or via mc++ <prefix>.in.xml which is preferable on a cluster.
