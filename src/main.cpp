#include <alps/parapack/parapack.h>
#include <alps/parapack/exchange.h>
int main(int argc,char** argv){
    try {
        return alps::parapack::start(argc,argv);
    }
    catch (std::exception& exc) {
        std::cerr << exc.what() << "\n";
            alps::comm_exit(true);
            return -1;
    }
    catch (...) {
        std::cerr << "Fatal Error: Unknown Exception!\n";
        return -2;
    }
}
