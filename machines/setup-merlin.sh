# RUN THIS SCRIPT ONLY THE FIRST TIME ON MERLIN!
# this clones the mc++ in the directory called, tries to install it and links the executable to ~/bin/mc++

# Load the required modules
source /opt/psi/config/profile.bash
module use unstable
module add gcc/4.7.4 openmpi/1.8.4 hdf5/1.8.14 boost/1.57.0
module add alps/2.2.b4

# clone the source
git clone https://github.com/domischi/mcpp.git

# install the file
cd mcpp/build
cmake ..
make

if [ -f ../bin/mc++ ]; then
    mkdir -p ~/bin
    ln -s ../bin/mc++ ~/bin/mc++
fi
