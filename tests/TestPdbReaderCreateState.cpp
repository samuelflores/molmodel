/* Regression test for PDBReader::createState fixes:
 *
 * 1. Multi-chain PDBs used to segfault: createState() indexed the
 *    compounds vector once per PDB chain, but createCompounds() only
 *    builds compounds for chains with recognizable biopolymer residues
 *    (waters/heterogens get no compound). A water chain after the
 *    protein chain therefore read past the end of the vector.
 *    It must instead skip compound-less chains with a warning.
 *
 * 2. Biopolymer::getResidue(name)/updResidue(name) used to dereference
 *    map::find(...)->second unchecked, crashing on unknown residue
 *    names. They must throw a descriptive exception instead.
 */

#include "molmodel/internal/Compound.h"
#include "molmodel/internal/PDBReader.h"
#include "molmodel/internal/Protein.h"

#include <iostream>
#include <fstream>
#include <cmath>

using namespace SimTK;
using namespace std;

// Minimal PDB: protein chain A (ALA-GLY) followed by a chain B of
// waters. Pre-fix, the water chain caused an out-of-bounds read of the
// compounds vector inside createState().
static const char* multichainPdb =
    "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
    "ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C\n"
    "ATOM      3  C   ALA A   1       1.980   1.430   0.000  1.00  0.00           C\n"
    "ATOM      4  O   ALA A   1       1.430   2.300   0.700  1.00  0.00           O\n"
    "ATOM      5  N   GLY A   2       3.180   1.560   0.100  1.00  0.00           N\n"
    "ATOM      6  CA  GLY A   2       3.850   2.840   0.150  1.00  0.00           C\n"
    "ATOM      7  C   GLY A   2       5.050   2.700   1.050  1.00  0.00           C\n"
    "ATOM      8  O   GLY A   2       5.900   3.600   1.100  1.00  0.00           O\n"
    "HETATM    9  O   HOH B   1      10.000  10.000  10.000  1.00  0.00           O\n"
    "HETATM   10  O   HOH B   2      12.000  10.000  10.000  1.00  0.00           O\n"
    "TER\n"
    "END\n";

void testCreateStateWithWaterChain() {
    // PDBReader's public API takes a filename only, so write the embedded
    // PDB to a temporary file.
    const std::string pdbFileName = "TestPdbReaderCreateState_tmp.pdb";
    {
        std::ofstream pdbFile(pdbFileName.c_str());
        assert(pdbFile.is_open());
        pdbFile << multichainPdb;
    }

    CompoundSystem system;
    SimbodyMatterSubsystem matter(system);
    DuMMForceFieldSubsystem forceField(system);
    forceField.loadAmber99Parameters();

    PDBReader pdb(pdbFileName, true);
    pdb.createCompounds(system, "");

    // The file has 2 chains (chain A protein + chain B waters), but only
    // the protein should have produced a compound.
    assert(system.getNumCompounds() == 1);

    system.modelCompounds();
    system.realizeTopology();
    State state = system.getDefaultState();

    // Pre-fix this segfaulted when it reached the water chain; it must
    // now skip the compound-less chain and return a fit quality.
    const Real rmsFit = pdb.createState(system, state);
    cout << "createState fit rms=" << rmsFit << endl;
    assert(std::isfinite(rmsFit));

    std::remove(pdbFileName.c_str());
}

void testUnknownResidueNameThrows() {
    // Pre-fix this dereferenced map::end() and crashed.
    Protein protein("AG");

    bool threw = false;
    try {
        protein.getResidue("9999");
    } catch (const std::exception& e) {
        threw = true;
        cout << "caught expected exception: " << e.what() << endl;
    }
    assert(threw);
}

int main() {
    try {
        testCreateStateWithWaterChain();
        testUnknownResidueNameThrows();
    }
    catch (const std::exception& e) {
        printf("EXCEPTION THROWN: %s\n", e.what());
        return 1;
    }
    catch (...) {
        printf("UNKNOWN EXCEPTION THROWN\n");
        return 1;
    }

    cout << "Done" << endl;
    return 0;
}
