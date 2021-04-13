/* vim: set sw=4 ts=4 sts=4 expandtab: */
/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Molmodel                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "molmodel/internal/PDBReader.h"
#include "molmodel/internal/Protein.h"
#include "molmodel/internal/RNA.h"
#include "molmodel/internal/DNA.h"
#include "SimTKsimbody.h"
#include <cctype>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <array>
#include <gemmi/cif.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/pdb.hpp>

#include "Representation.h"
#include "Util.h"

using namespace SimTK;
using std::map;
using std::string;
using std::vector;

static const std::array<std::string, 3> INPUT_FILE_EXTNS = { ".pdb", ".cif", ".cif.gz" };

static std::string getFileExtension( const std::string &filename ) {
    static const std::string DOT = ".";

    auto lastDot                 = filename.find_last_of( DOT );
    if ( lastDot == std::string::npos )
        return "";
    auto tail                    = filename.substr( lastDot );
    if ( tail == ".gz" ) {
        lastDot                  = filename.substr( 0, lastDot - 1 ).find_last_of( DOT );
        if (lastDot == std::string::npos )
            return ".gz";
        tail                     = filename.substr( lastDot );
    }

    std::transform               ( tail.begin(), tail.end(), tail.begin(), ::tolower );
    return tail;
}

static bool isCif( const std::string extension ) {
    return extension.find("cif") != std::string::npos;
}

static gemmi::Structure getStructureFromFile( const std::string &filename ) {
    if ( isCif ( getFileExtension( filename ) ) ) {
        //======================================== Read in the file into Gemmi document
        gemmi::cif::Document gemmiDoc                 = gemmi::cif::read ( gemmi::MaybeGzipped ( filename ) );
        //======================================== Check the number of blocks
        if ( gemmiDoc.blocks.size() < 1 )
        {
            std::cerr << "!!! Error !!! The file " << filename << " contains 0 blocks. Nothing the be read." << std::endl;
            std::exit                                 ( EXIT_FAILURE );
        }
        else if ( gemmiDoc.blocks.size() > 1 )
        {
            std::cerr << "!!! Warning !!! The file " << filename << " contains multiple blocks. Molmodel will use the first block named " << gemmiDoc.blocks.at(0).name << " and ignore the rest." << std::endl;
        }

        //======================================== Generate structure from block
        return gemmi::impl::make_structure_from_block ( gemmiDoc.blocks.at(0) );
    }

    // Assume PDB
    return gemmi::read_pdb_file                       ( filename );
}

class PDBReader::PDBReaderImpl {
public:
    PDBReaderImpl(string filename ) : hasBuiltSystem(false) { // second parameter is a vector of strings specifying chain ID, residue combinations to be deleted.  Optional parameter, defaults to an empty vector.
        std::cout<<__FILE__<<":"<<__LINE__<<"  filename.c_str()  >"<< filename.c_str()<<"< "<<std::endl;
        //============================================ Read in PDB or CIF
        const auto extension = getFileExtension(filename);

        if ( std::find_if( INPUT_FILE_EXTNS.cbegin(), INPUT_FILE_EXTNS.cend(), [&extension](const std::string &item) { return extension == item; } ) == INPUT_FILE_EXTNS.cend() ) {
            std::cerr << "!!! Error !!! This file type is not supported." << std::endl;
            std::exit ( EXIT_FAILURE );
        }

        //======================================== Generate structure from block
        gemmi::Structure gemmiStruct              = getStructureFromFile ( filename );

        //======================================== For each model
        int modNumGemmi = 0;
        for ( const auto &mod : gemmiStruct.models )
        {
            //==================================== Initialise the molmodel model
            std::string modName = "mol" + std::to_string(modNumGemmi++);

            Repr::Model rModel{std::move(modName), {}};

            //==================================== For each chain
            for ( const auto &chain : mod.chains)
            {
                Repr::Chain rChain{trim_both(chain.name)};

                //================================ For each residue
                int resNumGemmi                   = 0;
                for ( const auto &residue : chain.residues )
                {
                    Repr::Residue rResidue;

                    //============================ Get residue insertion code, residueID and residue name
                    char ICode                    = residue.seqid.icode;

                    if ( residue.seqid.num.has_value() ) { resNumGemmi = residue.seqid.num.value; }
                    else                                 { resNumGemmi++; }

                    rResidue.id = resNumGemmi;
                    rResidue.type = Repr::getResidueType(trim_both(residue.name));
                    rResidue.prop = Repr::getResidueProp(rResidue.type);
                    std::cout<<__FILE__<<":"<<__LINE__<<"res num = "<<rResidue.id<<" res type = "<<rResidue.type<< "residue name untrimmed = "<<residue.name<<"trimmmed residue name = >"<<trim_both(residue.name) <<"<"<<std::endl;
                    rResidue.insertion_code = ICode;

                    //============================ For each atom
                    for ( const auto &atom : residue.atoms )
                    {
                        //======================== Initialise atom related variables
                        Repr::Atom rAtom;

                        //======================== Get atom info
                        char altLoc               = atom.altloc;

                        //======================== Fill in atom information
                        rAtom.orig_id             = atom.serial;
                        rAtom.name                = Repr::AtomName{trim_both(atom.name)};
                        rAtom.resType             = rResidue.type;
                        rAtom.resProp             = rResidue.prop;

                        if ( !( altLoc == '\0' ) && !( altLoc == 'A' ) )
                            continue;                                                    // This is in keeping with the mol_DbPdbAtomProcLongChainId() function; it allows only the first alt-loc to be added.

                        rAtom.res_seq             = resNumGemmi;
                        rAtom.insertion_code      = ICode;
                        rAtom.posX                = atom.pos.x;
                        rAtom.posY                = atom.pos.y;
                        rAtom.posZ                = atom.pos.z;
                        rAtom.temp                = atom.b_iso;

                        rAtom.het                  = false;                       // We do not actually know, but let's go with false here ...

                        rAtom.setChainId(chain.name);

                        rResidue.atoms.emplace_back(std::move(rAtom));
                    }

                    rChain.residues.emplace_back(std::move(rResidue));
                }

                rModel.chains.emplace_back(std::move(rChain));
            }

            rStructure.emplace_back(std::move(rModel));
        }

        // NOTE:
        // There used to be a "structure count check"
        // As far as I can tell this was actually supposed to be "model count check"
        // Let us do the model count check
        std::cout<<__FILE__<<":"<<__LINE__<<" rStructure.size()  = "<< rStructure.size() <<std::endl;
        if (rStructure.size() > 1)
            throw new std::runtime_error{"Input file must contain only one model"};
        //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
    }

    void createCompounds( CompoundSystem& system, const String & chainsPrefix  ) {
        SimTK_APIARGCHECK_ALWAYS(!hasBuiltSystem, "PDBReaderImpl", "createSystem", "createSystem() has already been called");

        // Loop over chains and create a Biopolymer from each one.
        if (rStructure.size() > 0) {
            for (const auto &rChain : rStructure.front().chains) {
                // Create a string of the sequence.
                auto firstValidResidueType = Repr::ResidueType::UNKNOWN;
                string sequence;
                for (const auto &rResidue : rChain.residues) {
                    // Ignore anything that is not a canonal RNA or DNA residue
                    const auto type = rResidue.type;
                    if (Repr::residueIsDNA(type) || Repr::residueIsProtein(type) || Repr::residueIsRNA(type)) {
                        sequence += Repr::getResidueSpecifier(rResidue.type).shortName;

                        if (firstValidResidueType == Repr::ResidueType::UNKNOWN){
                            firstValidResidueType = rResidue.type;}
                    }
                }

                if (sequence.length() < 1){
                    std::cout<<__FILE__<<":"<<__LINE__<<""<<std::endl;
                    continue;}

                std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::toupper);

                if (Repr::residueIsRNA(firstValidResidueType)) {
                    std::cout<<__FILE__<<":"<<__LINE__<<"Creating an RNA"<<std::endl;
                    // Create an RNA.
                    RNA rna(sequence);
                    //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
                    for (int i = 0; i < rna.getNumResidues(); i ++) {
                        const auto &r = rChain.residues.at(i);
                        rna.updResidue(ResidueInfo::Index(i)).setPdbResidueNumber(r.id);
                        rna.updResidue(ResidueInfo::Index(i)).setPdbInsertionCode(r.insertion_code);//    ('101A');
                    }
                    rna.assignBiotypes();
                    std::cout<<__FILE__<<":"<<__LINE__<<" (*chains).id >"<<rChain.getId()<<"< "<<std::endl;
                    std::cout<<__FILE__<<":"<<__LINE__<<" created an RNA"<<std::endl;
                    std::cout<<__FILE__<<":"<<__LINE__<<" with chain >"<<rna.getPdbChainId()<<"< "<<std::endl;
                    compounds.push_back(rna);
                } else if (Repr::residueIsDNA(firstValidResidueType)) {
                    std::cout<<__FILE__<<":"<<__LINE__<<"Creating a DNA"<<std::endl;
                    // Create a  DNA.
                    DNA dna(sequence);
                    for (int i = 0; i < dna.getNumResidues(); i ++) {
                        const auto &r = rChain.residues.at(i);
                        dna.updResidue(ResidueInfo::Index(i)).setPdbResidueNumber(r.id);
                        dna.updResidue(ResidueInfo::Index(i)).setPdbInsertionCode(r.insertion_code);
                    }
                    dna.assignBiotypes();
                    compounds.push_back(dna);
                } else if (Repr::residueIsProtein(firstValidResidueType)) {
                    std::cout<<__FILE__<<":"<<__LINE__<<"Creating a protein"<<std::endl;
                    // Create a Protein.
                    // scf changed this to use one more parameter in the Protein constructor, set to "Torsion".  Default is "Rigid".  Now the peptide bond will not be rigid.
                    std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
                    Protein protein(sequence,BondMobility::Torsion);
                    //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;

                        protein.updResidue(ResidueInfo::Index(0)).setPdbResidueNumber(rChain.residues.at(0).id);
                        protein.updResidue(ResidueInfo::Index(0)).setPdbInsertionCode(' ');
                        //std::cout<<__FILE__<<":"<<__LINE__<<" i, residue number, insertion code, residue type : "<< protein.updResidue(ResidueInfo::Index(0)).getPdbResidueNumber()  <<", "<< protein.updResidue(ResidueInfo::Index(0)).getPdbInsertionCode()<<", "<<protein.updResidue(ResidueInfo::Index(0)).getOneLetterCode()<<std::endl;
                        //std::cout<<__FILE__<<":"<<__LINE__<<" i, residue number, insertion code, residue type : "<<chainResidues[0]->id<<", "<<chainResidues[0]->insertion_code<<", "<< protein.updResidue(ResidueInfo::Index(0)).getOneLetterCode()<<std::endl;
                    for (int i = 1; i < protein.getNumResidues() - 1; i++) { // assume protein capping is ON.  this means first and last residues are end caps, we ignore at this stage.
                        const auto &r = rChain.residues.at(i - 1);
                        protein.updResidue(ResidueInfo::Index(i)).setPdbResidueNumber(r.id);
                        protein.updResidue(ResidueInfo::Index(i)).setPdbInsertionCode(r.insertion_code);
                        //std::cout<<__FILE__<<":"<<__LINE__<<" i, residue number, insertion code, residue type : "<< protein.updResidue(ResidueInfo::Index(i)).getPdbResidueNumber()  <<", "<< protein.updResidue(ResidueInfo::Index(i)).getPdbInsertionCode()<<", "<<protein.updResidue(ResidueInfo::Index(i)).getOneLetterCode()<<std::endl;
                    }
                    // was unable to retrieve C-terminal end cap for some reason:
                        // protein.updResidue(ResidueInfo::Index((protein.getNumResidues() - 1) )).setPdbResidueNumber(chainResidues[(protein.getNumResidues() - 1) ]->id+1);
                        // protein.updResidue(ResidueInfo::Index((protein.getNumResidues() - 1) )).setPdbInsertionCode(' ');
                        //std::cout<<__FILE__<<":"<<__LINE__<<" i, residue number, insertion code, residue type : "<< protein.updResidue(ResidueInfo::Index((protein.getNumResidues() - 1) )).getPdbResidueNumber()  <<", "<< protein.updResidue(ResidueInfo::Index( (protein.getNumResidues() - 1) )).getPdbInsertionCode()<<", "<<protein.updResidue(ResidueInfo::Index( (protein.getNumResidues() - 1) )).getOneLetterCode()<<std::endl;
                    protein.assignBiotypes();
                    std::cout<<__FILE__<<":"<<__LINE__<<" (*chains).id >"<<rChain.getId()<<"< "<<std::endl;
                    std::cout<<__FILE__<<":"<<__LINE__<<" created an RNA"<<std::endl;
                    std::cout<<__FILE__<<":"<<__LINE__<<" with chain >"<<protein.getPdbChainId()<<"< "<<std::endl;
                    compounds.push_back(protein);
                } else {
                    const auto &r = rChain.residues.at(0);
                    std::cout<<__FILE__<<":"<<__LINE__<<" Did not recognize chainResidues[0]->type "<<r.type<<" ("<<Repr::getResidueSpecifier(r.type).longName<<" - "<<"). Please use only canonical RNA, DNA, and protein residue names"<<std::endl;
                    exit(EXIT_FAILURE);
                }
                std::cout<<__FILE__<<":"<<__LINE__<<" setPdbChainId(String("<<rChain.getId()<<") "<<std::endl;
                if (std::strlen(rChain.getId()) > 1)
                    std::cout << "!!! WARNING: Chain Id is longer than one character !!!" << std::endl;

                compounds[compounds.size()-1].setPdbChainId(rChain.getId());

                // Prepend chainId prefix if needed
                String chainId = rChain.getId();
                if ( !chainsPrefix.empty() ) {
                    chainId.insert( 0, chainsPrefix );
                    std::cout << __FILE__<< ":" << __LINE__ << " Chain: " << chainId << std::endl;
                    compounds[compounds.size() - 1].setPdbChainId(chainId);
                }
            }
	    }

        // Add them to the system.
        for (size_t i = 0; i < compounds.size(); i++)
            system.adoptCompound(compounds[i]);
        hasBuiltSystem = true;
    }

    Real createState(const CompoundSystem& system, State& state) const {
        SimTK_APIARGCHECK_ALWAYS(hasBuiltSystem, "PDBReaderImpl", "createState", "createSystem() has not yet been called");

        // Loop over atoms, match each one to the appropriate mobilized body, and
        // create a list of stations that will be used for fitting the State.

        map<MobilizedBodyIndex, vector<Vec3>> stations;
        map<MobilizedBodyIndex, vector<Vec3>> targetLocations;

        if (rStructure.size() > 0) {
            int proteinIndex = 0;
            for (const auto &rChain : rStructure.front().chains) {
                for (const auto &rResidue : rChain.residues) {
                    const auto &resName = std::to_string(rResidue.id);
                    const Biopolymer& compound = compounds[proteinIndex];
                    const ResidueInfo::Index resIx = compound.getResidue(resName).getIndex();
                    const ResidueInfo& residue = compound.getResidue(resName);

                    for (const auto &rAtom : rResidue.atoms) {
                        for (ResidueInfo::AtomIndex atomId(0); atomId < residue.getNumAtoms(); atomId++) {
                            if (rAtom.name.name == residue.getAtomName(atomId)) {
                                MobilizedBodyIndex bodyId = compound.getResidueAtomMobilizedBodyIndex(resIx, atomId);
                                stations[bodyId].push_back(compound.getResidueAtomLocationInMobilizedBodyFrame(resIx, atomId));
                                targetLocations[bodyId].push_back(Vec3(rAtom.posX, rAtom.posY, rAtom.posZ));
                            }
                        }
                    }
                }
                proteinIndex++;
            }
        }

        // Now perform the fitting.
        vector<MobilizedBodyIndex> bodyList;
        vector<vector<Vec3> > stationList;
        vector<vector<Vec3> > targetList;
        for (const auto &it : stations) {
            bodyList.push_back(it.first);
            stationList.push_back(it.second);
            targetList.push_back(targetLocations.find(it.first)->second);
        }
        // sherm 100307: Optimizers now use relative tolerance.
        Real tolerance = .001; // 0.1%
        return ObservedPointFitter::findBestFit(system, state, bodyList, stationList, targetList, tolerance);
    }

private:
    Repr::Structure rStructure;
    vector<Biopolymer> compounds;
    bool hasBuiltSystem;
};

PDBReader::PDBReader(string filename ) {
    //std::cout<<__FILE__<<":"<<__LINE__<<" >"<< deletedResidueVector.size() <<"<"<<std::endl;
    impl = new PDBReaderImpl(filename );
}

PDBReader::~PDBReader() {
    delete impl;
}

void PDBReader::createCompounds(CompoundSystem& system, const String & chainsPrefix = "" ) {
    impl->createCompounds(system, chainsPrefix);
}

Real PDBReader::createState(const CompoundSystem& system, State& state) const {
    return impl->createState(system, state);
}

