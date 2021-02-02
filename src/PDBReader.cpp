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
#include "mol.h"
#include "SimTKsimbody.h"
#include <cctype>
#include <map>
#include <vector>
#include <algorithm>

#ifdef GEMMI_USAGE
    //================================================ Include Gemmi headers
    #include <gemmi/gz.hpp>
#endif

using namespace SimTK;
using std::map;
using std::string;
using std::vector;

class PDBReader::PDBReaderImpl {
public:
    PDBReaderImpl(string filename ) : hasBuiltSystem(false) { // second parameter is a vector of strings specifying chain ID, residue combinations to be deleted.  Optional parameter, defaults to an empty vector.
        
        std::cout<<__FILE__<<":"<<__LINE__<<"  filename.c_str()  >"<< filename.c_str()<<"< "<<std::endl;
        
        //============================================ Read in PDB or CIF
        if ( filename.substr ( filename.length() - 4, filename.length() - 1) == ".pdb" )
        {
            mol_DbRead("mol", filename.c_str(), MOL_DB_PDB, &model );
        }
        else if ( ( filename.substr ( filename.length() - 4, filename.length() - 1) == ".cif" ) ||
                  ( filename.substr ( filename.length() - 7, filename.length() - 1) == ".cif.gz" ) )
        {
#ifdef GEMMI_USAGE
            //======================================== Read in the file into Gemmi document
            gemmi::cif::Document gemmiDoc             = gemmi::cif::read ( gemmi::MaybeGzipped ( filename ) );
            
            //======================================== Check the number of blocks
            if ( gemmiDoc.blocks.size() < 1 )
            {
                std::cerr << "!!! Error !!! The file " << filename << " contains 0 blocks. Nothing the be read." << std::endl;
                exit                                  ( EXIT_FAILURE );
            }
            else if ( gemmiDoc.blocks.size() > 1 )
            {
                std::cerr << "!!! Warning !!! The file " << filename << " contains multiple blocks. Molmodel will use the first block named " << gemmiDoc.blocks.at(0).name << " and ignore the rest." << std::endl;
            }

            //======================================== Generate structure from block
            gemmi::Structure gemmiStruct              = gemmi::impl::make_structure_from_block ( gemmiDoc.blocks.at(0) );
            
            //======================================== For each model
            for ( unsigned int moIt = 0; moIt < static_cast<unsigned int> ( gemmiStruct.models.size() ); moIt++ )
            {
                //==================================== Initialise the molmodel model
                std::stringstream molName;
                molName << "mol" << moIt;
                mol_MolModelCreate                    ( molName.str().c_str(), &model );
                
                //==================================== For each chain
                for ( unsigned int chIt = 0; chIt < static_cast<unsigned int> ( gemmiStruct.models.at(moIt).chains.size() ); chIt++ )
                {
                    //================================ Get the chain ID
                    std::string chainId               =  std::string ( gemmiStruct.models.at(moIt).chains.at(chIt).name );
                 
                    //================================ For each residue
                    int resNumGemmi                   = 0;
                    for ( unsigned int reIt = 0; reIt < static_cast<unsigned int> ( gemmiStruct.models.at(moIt).chains.at(chIt).residues.size() ); reIt++ )
                    {
                        //============================ Get residue insertion code, residueID and residue name
                        char ICode                    = gemmiStruct.models.at(moIt).chains.at(chIt).residues.at(reIt).seqid.icode;
                        
                        if ( gemmiStruct.models.at(moIt).chains.at(chIt).residues.at(reIt).seqid.num.has_value() ) { resNumGemmi = gemmiStruct.models.at(moIt).chains.at(chIt).residues.at(reIt).seqid.num.value; }
                        else                                                                                       { resNumGemmi++; }
                        
                        PdbResidueId residueId        ( resNumGemmi, ICode );
                        std::string residueName       = gemmiStruct.models.at(moIt).chains.at(chIt).residues.at(reIt).name;
                       
                        //============================ For each atom
                        for ( unsigned int atIt = 0; atIt < static_cast<unsigned int> ( gemmiStruct.models.at(moIt).chains.at(chIt).residues.at(reIt).atoms.size() ); atIt++ )
                        {
                            //======================== Initialise atom related variables
                            MolAtom atom;
                            MolStructure *struc;
                            
                            //======================== Get atom info
                            char altLoc               = gemmiStruct.models.at(moIt).chains.at(chIt).residues.at(reIt).atoms.at(atIt).altloc;
                            
                            //======================== Fill in atom information
                            atom.orig_id              = gemmiStruct.models.at(moIt).chains.at(chIt).residues.at(reIt).atoms.at(atIt).serial;
                            atom.name                 = mol_StrCopy( gemmiStruct.models.at(moIt).chains.at(chIt).residues.at(reIt).atoms.at(atIt).name.c_str(), model );
                            mol_ResTypeConv           ( residueName.c_str(), &atom.res_type );
                            atom.res_prop             = mol_res_props[atom.res_type];
                            
                            if ( !( altLoc == '\0' ) && !( altLoc == 'A' ) )
                                continue;                                                    // This is in keeping with the mol_DbPdbAtomProcLongChainId() function; it allows only the first alt-loc to be added.
                            
                            atom.res_seq              = resNumGemmi;
                            atom.insertion_code       = ICode;
                            atom.pos[0]               = gemmiStruct.models.at(moIt).chains.at(chIt).residues.at(reIt).atoms.at(atIt).pos.x;
                            atom.pos[1]               = gemmiStruct.models.at(moIt).chains.at(chIt).residues.at(reIt).atoms.at(atIt).pos.y;
                            atom.pos[2]               = gemmiStruct.models.at(moIt).chains.at(chIt).residues.at(reIt).atoms.at(atIt).pos.z;
                            atom.temp                 = gemmiStruct.models.at(moIt).chains.at(chIt).residues.at(reIt).atoms.at(atIt).b_iso;
                            
                            atom.het                  = MOL_FALSE;                       // We do not actually know, but let's go with false here ...
                            
                            if (atom.name[0] != ' ')  { mol_AtomTypeConv ( &atom.name[0], &atom.type ); }
                            else                      { mol_AtomTypeConv ( &atom.name[1], &atom.type ); }
                            
                            atom.chain_id             = chainId[0];
                            atom.long_chain_id        = mol_StrCopy( " ", model );       // We do not have LongChainID, so empty
                            
                            mol_MolModelCurrStrucGet  ( model, &struc );                 // Copies model->curr_struc to struc
                            mol_StructureAtomAdd      ( struc, false, &atom );           // Saves the atom to the model (through struc)
                        }
                    }
                }
                
                //==================================== Build biomolecule from the filled in structs
                MolStructure *struc;
                mol_MolModelCurrStrucGet              ( model, &struc ); // Fill in the struc structure
                mol_StructureChainsBuild              ( struc, 1 );      // Build protein chains
                mol_StructureChainsBuild              ( struc, 2 );      // Build solvent chains
            }
#else
            std::cerr << "ERROR: The suplied file has the mmCIF format, but Molmodel was not compiled with the Gemmi library, which is required for the mmCIF support. Please re-run the Molmodel cmake command with the -DUSE_GEMMI=TRUE -DGEMMI_PATH=/path/to/gemmi/include option (and then alse re-run the make install command). Terminating..." << std::endl;
            exit                                      ( -1 );
#endif
        }
        else
        {
            std::cerr << "Failed to detect the extension of the input file " << filename << ". The supported extensions are: \'.pdb\' and \'.cif\' (or \'.cif.gz\'). Terminating now..." << std::endl;
            exit                                      ( -1 );
        }

        //============================================ Proceed as normal
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        MolStructure *structure;
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        int numStructures;
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        mol_MolModelStructuresGet ( model, &numStructures, &structure );
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
        SimTK_APIARGCHECK_ALWAYS(numStructures == 1, "PDBReaderImpl", "PDBReaderImpl", "The PDB file must contain exactly one structure");
        std::cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<std::endl;
    }

    ~PDBReaderImpl() {
        mol_MemFreeModel(model);
    }

    void createCompounds( CompoundSystem& system, const String & chainsPrefix  ) {
        SimTK_APIARGCHECK_ALWAYS(!hasBuiltSystem, "PDBReaderImpl", "createSystem", "createSystem() has already been called");
        MolStructure *structure;
        int numStructures;
        mol_MolModelStructuresGet(model, &numStructures, &structure); // numStructures and structure are just copied from corresponding members of model.  But where does model come from?  Apparently it's a private member of PDBReader, set by mol_DbRead above
        
        //std::vector<MolAtom> atoms;
        MolAtom* atoms;
        int numAtoms;
        mol_StructureAtomsGet(structure, &numAtoms, &atoms);
        MolChain* chains;
        int numChains;
        mol_StructureChainsGet(structure, &numChains, &chains); // scf: members of structure are simply copied to numChains and chains.
        
        // Loop over chains and create a Biopolymer from each one.
        while (chains) {
            MolResidue** chainResidues;
            int numResidues;
            mol_ChainResiduesGet(chains, &numResidues, &chainResidues);
            
            // Create a string of the sequence.
            
            string sequence;
            for (int i = 0; i < numResidues; ++i) {
                std::cout<<__FILE__<<":"<<__LINE__<<" Adding a residue of type >"<< mol_res_names[chainResidues[i]->type][2]<<"<"  <<std::endl;
                sequence += mol_res_names[chainResidues[i]->type][2];
            }
            std::transform(sequence.begin(), sequence.end(), sequence.begin(), (int(*)(int)) std::toupper);
            if ((chainResidues[0]->type >= 25) &&
                (chainResidues[0]->type <= 28)) 
            {
                std::cout<<__FILE__<<":"<<__LINE__<<"Creating an RNA"<<std::endl;
                // Create an RNA.
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
                RNA rna(sequence);
                //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
                for (int i = 0; i < rna.getNumResidues(); i ++) {
                    rna.updResidue(ResidueInfo::Index(i)).setPdbResidueNumber(chainResidues[i]->id);
                    rna.updResidue(ResidueInfo::Index(i)).setPdbInsertionCode(chainResidues[i]->insertion_code);//    ('101A');
                }
                rna.assignBiotypes();
                std::cout<<__FILE__<<":"<<__LINE__<<" (*chains).name >"<<(*chains).name<<"< "<<std::endl;
                std::cout<<__FILE__<<":"<<__LINE__<<" (*chains).id >"<<(*chains).id<<"< "<<std::endl;
                std::cout<<__FILE__<<":"<<__LINE__<<" (*chains).long_chain_id >"<<(*chains).long_chain_id<<"< "<<std::endl;
                std::cout<<__FILE__<<":"<<__LINE__<<" created an RNA"<<std::endl; 
                std::cout<<__FILE__<<":"<<__LINE__<<" with chain >"<<rna.getPdbChainId()<<"< "<<std::endl;
                compounds.push_back(rna);
            } 
            else if ((chainResidues[0]->type >= 40) &&
                     (chainResidues[0]->type <= 43)) {
                std::cout<<__FILE__<<":"<<__LINE__<<"Creating a DNA"<<std::endl;
                // Create a  DNA.
                
                DNA dna(sequence);
                for (int i = 0; i < dna.getNumResidues(); i ++) {
                    dna.updResidue(ResidueInfo::Index(i)).setPdbResidueNumber(chainResidues[i]->id);
                    dna.updResidue(ResidueInfo::Index(i)).setPdbInsertionCode(chainResidues[i]->insertion_code);
                }
                dna.assignBiotypes();
                compounds.push_back(dna);
            }  
            else if (((chainResidues[0]->type < 21) && (chainResidues[0]->type > 0)) ||  // type 0 is "unknown".  1-20 are protein amino acids
                     (chainResidues[0]->type == 44) ) {                                  // type 44 is CYX
                std::cout<<__FILE__<<":"<<__LINE__<<"Creating a protein"<<std::endl;
                // Create a Protein.
                // scf changed this to use one more parameter in the Protein constructor, set to "Torsion".  Default is "Rigid".  Now the peptide bond will not be rigid.
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
                Protein protein(sequence,BondMobility::Torsion);
                //std::cout<<__FILE__<<":"<<__LINE__<<std::endl;

                    protein.updResidue(ResidueInfo::Index(0)).setPdbResidueNumber(chainResidues[0]->id-1);
                    protein.updResidue(ResidueInfo::Index(0)).setPdbInsertionCode(' ');
                    //std::cout<<__FILE__<<":"<<__LINE__<<" i, residue number, insertion code, residue type : "<< protein.updResidue(ResidueInfo::Index(0)).getPdbResidueNumber()  <<", "<< protein.updResidue(ResidueInfo::Index(0)).getPdbInsertionCode()<<", "<<protein.updResidue(ResidueInfo::Index(0)).getOneLetterCode()<<std::endl;
                    //std::cout<<__FILE__<<":"<<__LINE__<<" i, residue number, insertion code, residue type : "<<chainResidues[0]->id<<", "<<chainResidues[0]->insertion_code<<", "<< protein.updResidue(ResidueInfo::Index(0)).getOneLetterCode()<<std::endl;
                for (int i = 1; i < (protein.getNumResidues() - 1); i ++) { // assume protein capping is ON.  this means first and last residues are end caps, we ignore at this stage.
                    protein.updResidue(ResidueInfo::Index(i)).setPdbResidueNumber(chainResidues[i-1]->id);
                    protein.updResidue(ResidueInfo::Index(i)).setPdbInsertionCode(chainResidues[i-1]->insertion_code);
                    //std::cout<<__FILE__<<":"<<__LINE__<<" i, residue number, insertion code, residue type : "<< protein.updResidue(ResidueInfo::Index(i)).getPdbResidueNumber()  <<", "<< protein.updResidue(ResidueInfo::Index(i)).getPdbInsertionCode()<<", "<<protein.updResidue(ResidueInfo::Index(i)).getOneLetterCode()<<std::endl;
                }
	            // was unable to retrieve C-terminal end cap for some reason:	 	
                    // protein.updResidue(ResidueInfo::Index((protein.getNumResidues() - 1) )).setPdbResidueNumber(chainResidues[(protein.getNumResidues() - 1) ]->id+1);
                    // protein.updResidue(ResidueInfo::Index((protein.getNumResidues() - 1) )).setPdbInsertionCode(' ');
                    //std::cout<<__FILE__<<":"<<__LINE__<<" i, residue number, insertion code, residue type : "<< protein.updResidue(ResidueInfo::Index((protein.getNumResidues() - 1) )).getPdbResidueNumber()  <<", "<< protein.updResidue(ResidueInfo::Index( (protein.getNumResidues() - 1) )).getPdbInsertionCode()<<", "<<protein.updResidue(ResidueInfo::Index( (protein.getNumResidues() - 1) )).getOneLetterCode()<<std::endl;
                protein.assignBiotypes();
                std::cout<<__FILE__<<":"<<__LINE__<<" (*chains).name >"<<(*chains).name<<"< "<<std::endl;
                std::cout<<__FILE__<<":"<<__LINE__<<" (*chains).id >"<<(*chains).id<<"< "<<std::endl;
                std::cout<<__FILE__<<":"<<__LINE__<<" (*chains).long_chain_id >"<<(*chains).long_chain_id<<"< "<<std::endl;
                std::cout<<__FILE__<<":"<<__LINE__<<" created an RNA"<<std::endl; 
                std::cout<<__FILE__<<":"<<__LINE__<<" with chain >"<<protein.getPdbChainId()<<"< "<<std::endl;
                compounds.push_back(protein);
            } else {
                std::cout<<__FILE__<<":"<<__LINE__<<" Did not recognize chainResidues[0]->type "<<chainResidues[0]->type<<". Please use only canonical RNA, DNA, and protein residue names"<<std::endl;
                exit(1);
            }   
            std::cout<<__FILE__<<":"<<__LINE__<<" setPdbChainId(String("<<(*chains).long_chain_id <<") "<<std::endl;
            if (std::string((*chains).long_chain_id).compare(std::string(" ")) != 0) {// if we are using a long chain ID
                std::cout<<__FILE__<<":"<<__LINE__<<" (*chains).long_chain_id = >"<<(*chains).long_chain_id<<"< "<<std::endl;
                std::cout<<__FILE__<<":"<<__LINE__<<" it would appear  a long chain ID is desired"<<std::endl;
                std::cout<<__FILE__<<":"<<__LINE__<<" (*chains).id = >"<<(*chains).id<<"< "<<std::endl;
                if (((*chains).id) != ((' ')) ) {
                    std::cout<<__FILE__<<":"<<__LINE__<<" But your PDB chain ID is not \" \" as it should be! "<<std::endl; exit(1);
                }
            }
            if (std::string((*chains).long_chain_id).compare(std::string(" ")) != 0) {
                std::cout<<__FILE__<<":"<<__LINE__<<" Detected that your (*chains).long_chain_id = >"<<(*chains).long_chain_id<<"< is NOT a blank space, >"<<std::string(" ")<<"< "<<std::endl;
                compounds[compounds.size()-1].setPdbChainId(String((*chains).long_chain_id ));
                std::cout<<__FILE__<<":"<<__LINE__<<std::endl;//" Detected that your (*chains).long_chain_id = >"<<(*chains).long_chain_id<<"< is NOT a blank space, >"<<std::string(" ")<<"< "<<std::endl;
            } else {
                std::cout<<__FILE__<<":"<<__LINE__<<" Detected that your (*chains).long_chain_id = >"<<(*chains).long_chain_id<<"< is a blank space, >"<<std::string(" ")<<"< "<<  std::endl;
                compounds[compounds.size()-1].setPdbChainId(String((*chains).id ));
            }
            // Prepend chainId prefix if needed
            String chainId = String((*chains).id);
            if ( !chainsPrefix.empty() ) {
                chainId.insert( 0, chainsPrefix );
                std::cout << __FILE__<< ":" << __LINE__ << " Chain: " << chainId << std::endl;
                compounds[compounds.size()-1].setPdbChainId( chainId );
            }
            chains = chains->next;
        }
        
        // Add them to the system.
        
        for (int i = 0; i < (int)compounds.size(); ++i)
            system.adoptCompound(compounds[i]);
        hasBuiltSystem = true;
    }
    
    Real createState(const CompoundSystem& system, State& state) const {
        SimTK_APIARGCHECK_ALWAYS(hasBuiltSystem, "PDBReaderImpl", "createState", "createSystem() has not yet been called");

        MolStructure *structure;
        int numStructures;
        mol_MolModelStructuresGet(model, &numStructures, &structure);
        //std::vector<MolAtom> atoms;
        MolAtom* atoms;
        int numAtoms;
        mol_StructureAtomsGet(structure, &numAtoms, &atoms);
        MolChain* chains;
        int numChains;
        mol_StructureChainsGet(structure, &numChains, &chains);

        // Loop over atoms, match each one to the appropriate mobilized body, and
        // create a list of stations that will be used for fitting the State.

        map<MobilizedBodyIndex, vector<Vec3> > stations;
        map<MobilizedBodyIndex, vector<Vec3> > targetLocations;
        int proteinIndex = 0;
        while (chains) {
            MolResidue** chainResidues;
            int numResidues;
            mol_ChainResiduesGet(chains, &numResidues, &chainResidues);
            for (int res = 0; res < numResidues; ++res) {
                MolAtomList resAtoms = chainResidues[res]->atoms;
                char resName[32];
                sprintf(resName, "%d", res);
                const Biopolymer& compound = compounds[proteinIndex];
                const ResidueInfo::Index resIx = compound.getResidue(resName).getIndex();
                const ResidueInfo& residue = compound.getResidue(resName);
// residue.setBondMobility(
                for (int i = 0; i < resAtoms.num; ++i) {
                    MolAtom atom = atoms[resAtoms.list[i]];
                    string atomName = atom.name;
                    atomName = atomName.erase(atomName.find_last_not_of(' ')+1);
                    atomName.erase(0, atomName.find_first_not_of(' '));
                    for (ResidueInfo::AtomIndex atomId(0); atomId < residue.getNumAtoms(); ++atomId) {
                        if (atomName == residue.getAtomName(atomId)) {
                            MobilizedBodyIndex bodyId = compound.getResidueAtomMobilizedBodyIndex(resIx, atomId);
                            stations[bodyId].push_back(compound.getResidueAtomLocationInMobilizedBodyFrame(resIx, atomId));
                            targetLocations[bodyId].push_back(Vec3(atom.pos[0], atom.pos[1], atom.pos[2]));
                        }
                    }
                }
            }
            proteinIndex++;
            chains = chains->next;
        }

        // Now perform the fitting.
        
        vector<MobilizedBodyIndex> bodyList;
        vector<vector<Vec3> > stationList;
        vector<vector<Vec3> > targetList;
        for (map<MobilizedBodyIndex, vector<Vec3> >::const_iterator iter = stations.begin(); iter != stations.end(); iter++) {
            bodyList.push_back(iter->first);
            stationList.push_back(iter->second);
            targetList.push_back(targetLocations.find(iter->first)->second);
        }
        // sherm 100307: Optimizers now use relative tolerance.
        Real tolerance = .001; // 0.1%
        return ObservedPointFitter::findBestFit(system, state, bodyList, stationList, targetList, tolerance);
    }
        
private:
    MolModel *model;
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

