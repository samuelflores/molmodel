#include "SimTKmolmodel.h"
#include <exception>
#include <iostream>
#include <fstream>

using namespace SimTK;
using namespace std;

#define ASSERT(x, msg) SimTK_ASSERT_ALWAYS(x, msg)

int main(int argc, char **argv) {
    ASSERT(argc == 2, "Invalid number of arguments, expected 2");

    ifstream tinkerStream(argv[1]);
    ASSERT(tinkerStream.good(), "Cannot read TinkerAmber99 parameters file");

    try {
        Biotype argonBiotype = Biotype::Argon();
        DuMMForceFieldSubsystem dumm;
        dumm.populateFromTinkerParameterFile(tinkerStream);

        Biotype::generateAllBiotypeCode(cout);
        return 0;
    } catch (const exception &ex) {
        cerr << "EXCEPTION THROWN: " << ex.what() << endl;
        return 1;
    } catch (...) {
        cerr << "UNKNOWN EXCEPTION THROWN" << endl;
        return 1;
    }
}
