This is mc++, a program mainly used to simulate XY systems in 2D with dipolar interaction. This is the basis for the PhD of Dominik Schildknecht <dominik.schildknecht@psi.ch>
The program is not part of the ALPS project (http://alps.comp-phys.org/), however is designed to be linked against it. However for anything surrounding the main simulation as input file generation and evaluation of the raw data alpspython is used. The software is therefore under the ALPS licencse v1.1, distributed with this program.

The papers needed to cite to use this program are displayed in the beginning of a run. Please note that it is mandatory to cite these if the results are published.

Installation Guide:
    - Download the ALPS library (http://alps.comp-phys.org/) and extract it
    - Set the environment variable ALPS_DIR to where ever you installed the library to (e.g. /opt/alps/share/alps)
    - Clone the mc++ program
        git clone TODO get path
    - Change into mc++ folder
        cd mc++
    - Go into the build directory and use cmake & make
        cd build; cmake ..; make

User Guide:
    Follow mainly the example. If so you are ensured to generate files that are compatible with this program. You can either start the program directly in python as given in the examples, or via mc++ <prefix>.in.xml which is preferable on a cluster.
