#include "SimTKmolmodel.h"
#include <iostream>
#include <string>
#include <stdexcept>

using namespace std;
using namespace SimTK;

void testCompletePdbLineShouldBeOK() 
{
    // Truncate at character 57; has coordinates but nothing after
    istringstream completeLine("ATOM    866  CA  LEU A 112      30.133  37.313  21.927  0.90 34.39           C  ");
    PdbStructure pdbStructure(completeLine, PdbStructure::InputType::PDB);

    const PdbAtom& atom = pdbStructure.getAtom(" CA ", PdbResidueId(112, ' '), "A");
    SimTK_TEST_EQ_TOL(atom.getOccupancy(), 0.90, 1e-7);
    SimTK_TEST_EQ_TOL(atom.getTemperatureFactor(),  34.39, 1e-7);
}

// Lack of coordinates should raise exception
void testMissingZCoordinateShouldRaiseException() 
{
    try {
        istringstream trunc("ATOM      1  N   ALA A   1     -52.630  -1.437");
        PdbStructure pdbStructure(trunc, PdbStructure::InputType::PDB);
    }
    catch (...) {
        return;
    }

    throw std::logic_error("Missing z-coordinate failed to raise an exception");
}

int main() 
{
try {
    testCompletePdbLineShouldBeOK();
    testMissingZCoordinateShouldRaiseException();

    cout << "PASSED" << endl;
    return 0;
}

catch (const std::exception& e)
{
    cerr << "EXCEPTION THROWN: " << e.what() << endl;

    cerr << "FAILED" << endl;
    return 1;
}

catch (...)
{
    cerr << "UNKNOWN EXCEPTION THROWN" << endl;

    cerr << "FAILED" << endl;
    return 1;
}

}

