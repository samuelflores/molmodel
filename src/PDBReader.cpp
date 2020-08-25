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

#ifdef MMDB2_LIB_USAGE
  #include <mmdb2/mmdb_manager.h>
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
#ifdef MMDB2_LIB_USAGE
            //======================================== Open file
            mmdb::CoorManager *mfile                  = new mmdb::CoorManager ( );
            
            //======================================== Check that the file was opened correctly
            if ( mfile->ReadCoorFile ( filename.c_str() ) )
            {
                std::cout << "MMDB2 Failed to open the file " << filename << ". Terminating now..." << std::endl;
                exit                                  ( -1 );
            }
        
            //======================================== Initialise the molmodel model
            mol_MolModelCreate                        ( "mol", &model );
            
            //======================================== Initialise MMDB crawl
            int noModels                              = 0;
            int noChains                              = 0;
            int noResidues                            = 0;
            int noAtoms                               = 0;
            mmdb::Model **mmdb2Model;
            mmdb::Chain **mmdb2Chain;
            mmdb::Residue **mmdb2Residue;
            mmdb::Atom **mmdb2Atom;
            
            //======================================== Find all models in the PDB file
            mfile->GetModelTable                      ( mmdb2Model, noModels );
            
            //======================================== Check that at least one model was found
            if ( noModels < 1 )
            {
                std::cerr << "MMDB2 found no models the file " << filename << ". Terminating now..." << std::endl;
                exit                                  ( -1 );
            }
            
            //======================================== If there are multiple models, just the first one will be used. This may needs some tweaking in the future
            if ( noModels < 1 )
            {
                std::cout << "WARNING: MMDB2 found multiple models (" << noModels << ") the file " << filename << ". Only the first model will be used." << std::endl;
            }
            
            //======================================== Check the first model is readable
            if ( mmdb2Model[0] )
            {
                //==================================== Deal with secondary structure's (alpha-helices first)
                int noHelices                         = mmdb2Model[0]->GetNumberOfHelices ( );
                
                for ( int hNo = 0; hNo < noHelices; hNo++ )
                {
                    //================================ Initialise variables
                    MolStructure *hStruc;
                    MolHelix helix;
                    mmdb::Helix *mmdb2Helix           = mmdb2Model[0]->GetHelix ( hNo+1 );
                    
                    //================================ Copy relevant information to molmodel structures
                    helix.num                         = mmdb2Helix->serNum;
                    mol_DbRecordStrGet                ( mmdb2Helix->helixID,     1, 3, helix.id );
                    mol_DbRecordStrGet                ( mmdb2Helix->initResName, 1, 3, helix.init_res_name );
                    helix.init_chain_id               = mmdb2Helix->initChainID[0];
                    helix.init_seq_num                = mmdb2Helix->initSeqNum;
                    mol_DbRecordStrGet                ( mmdb2Helix->endResName,  1, 3, helix.term_res_name );
                    helix.term_chain_id               = mmdb2Helix->endChainID[0];
                    helix.term_seq_num                = mmdb2Helix->endSeqNum;
                    helix.helixClass                  = mmdb2Helix->helixClass;
                    helix.length                      = mmdb2Helix->length;
                    
                    //================================ Save the read helix to the molmodel model
                    mol_MolModelCurrStrucGet          ( model, &hStruc ); // Copies model->curr_struc to hStruc
                    
                    MolHelix *hlist;
                    int num                           = hStruc->secondary.helix.num;
                    int size                          = hStruc->secondary.helix.size;
                    hlist                             = hStruc->secondary.helix.list;
                    if (num == (size - 1))
                    {
                        size                         += 1000;
                        mem_Realloc                   ( hlist, size, hStruc->model, MolHelix* );
                        hStruc->secondary.helix.size  = size;
                        hStruc->secondary.helix.list  = hlist;
                    }
                    memcpy                            ( hlist+num, &helix, sizeof(MolHelix) );
                    num                              += 1;
                    hStruc->secondary.helix.num       = num;
                }
                
                //==================================== Deal with secondary structure's (beta-sheets now)
                int noSheets                          = mmdb2Model[0]->GetNumberOfSheets ( );
                
                for ( int sNo = 0; sNo < noSheets; sNo++ )
                {
                    //================================ Initialise variables
                    mmdb::Sheet *mmdb2Sheet           = mmdb2Model[0]->GetSheet ( sNo+1 );
                    
                    int noStrands                     = mmdb2Sheet->nStrands;
                    for ( int stNo = 0; stNo < noStrands; stNo++ )
                    {
                        //============================ Initialise variables
                        MolStructure *sStruct;
                        MolSheet sheet;
                        mmdb::Strand *mmdb2Strand     = mmdb2Sheet->strand[stNo];
                        
                        //============================ Copy relevant information to molmodel structures
                        sheet.num                     = mmdb2Strand->strandNo;
                        mol_DbRecordStrGet            ( mmdb2Strand->sheetID,     1, 2, sheet.id );
                        sheet.num_strands             = noStrands;
                        mol_DbRecordStrGet            ( mmdb2Strand->initResName, 1, 3, sheet.init_res_name ); 
                        sheet.init_chain_id           = mmdb2Strand->initChainID[0];
                        sheet.init_seq_num            = mmdb2Strand->initSeqNum;
                        sheet.init_icode              = mmdb2Strand->initICode[0];
                        mol_DbRecordStrGet            ( mmdb2Strand->endResName,  1, 3, sheet.term_res_name );
                        sheet.term_chain_id           = mmdb2Strand->endChainID[0];
                        sheet.term_seq_num            = mmdb2Strand->endSeqNum;
                        sheet.end_icode               = mmdb2Strand->endICode[0];
                        sheet.sense                   = mmdb2Strand->sense;
                        
                        //============================ Save the read sheet to the molmodel model
                        mol_MolModelCurrStrucGet      ( model, &sStruct ); // Copies model->curr_struc to sStructt
                        
                        MolStructureSecSheet *secStruct;
                        MolSheet *list;
                        secStruct                     = &sStruct->secondary.sheet;
                        
                        int num                       = secStruct->num;
                        int size                      = secStruct->size;
                        list                          = secStruct->list;
                        if ( num == (size - 1) )
                        {
                          size                       += 1000;
                          mem_Realloc                 ( list, size, sStruct->model, MolSheet* );
                          secStruct->size             = size;
                          secStruct->list             = list;
                        }

                        memcpy                        ( list+num, &sheet, sizeof(MolSheet) );
                        num                          += 1;
                        secStruct->num                = num;
                    }
                    
                }
                
                //==================================== Find all chains for first model (the 1 in the first argument to the GetCHainTable() function)
                mfile->GetChainTable                  ( 1, mmdb2Chain, noChains );
                
                //==================================== For each chain
                for ( int nCh = 0; nCh < noChains; nCh++ )
                {
                    //================================ Check that the chain is reaable
                    if ( mmdb2Chain[nCh] )
                    {
                        //============================ Get all residues for this chain
                        mfile->GetResidueTable        ( 1, nCh, mmdb2Residue, noResidues );
                        
                        //============================ For each residue
                        for ( int nRes = 0; nRes < noResidues; nRes++ )
                        {
                            //======================== Check that the residue is readable
                            if ( mmdb2Residue[nRes] )
                            {
                                //==================== Get all atoms
                                mfile->GetAtomTable   ( 1, nCh, nRes, mmdb2Atom, noAtoms );
                                
                                //==================== For each atom
                                for ( int aNo = 0; aNo < noAtoms; aNo++ )
                                {
                                    //================ Check that the atom is readable
                                    if ( mmdb2Atom[aNo] )
                                    {
                                        //============ Check for termination 'residue'
                                        if ( mmdb2Atom[aNo]->Ter )
                                        {
                                            //======== Initialise variables
                                            MolStructure *terStruc;
                                            MolChainTerm term;
                                            
                                            term.id   = mmdb2Atom[aNo]->serNum;
                                            mol_DbRecordStrGet ( mmdb2Residue[nRes]->name, 1, 3, term.res_name );
                                            term.chain_id = mmdb2Chain[nCh]->GetChainID()[0];
                                            term.res_seq = mmdb2Residue[nRes]->seqNum;
                                            std::string ICodeHlp = std::string ( mmdb2Residue[nRes]->GetInsCode() );
                                            char ICode;
                                            if ( ICodeHlp == "" ) { ICode = ' '; } else { ICode = mmdb2Residue[nRes]->GetInsCode()[0]; }
                                            term.insertion_code = ICode;
                                            
                                            //============================ Save the read TER to the molmodel model
                                            mol_MolModelCurrStrucGet ( model, &terStruc );
                                            
                                            MolChainTerm *tp;
                                            mem_Alloc ( tp, 1, model, MolChainTerm* );
                                            tp->id    = term.id;
                                            tp->chain_id = term.chain_id;
                                            tp->res_seq = term.res_seq;
                                            tp->insertion_code = term.insertion_code;
                                            strcpy    ( tp->res_name, term.res_name );
                                            tp->next = terStruc->terms;
                                            terStruc->terms = tp;
                                            
                                            continue;
                                        }
                                        
                                        //============ Process atom
                                        MolAtom atom;
                                        MolStructure *struc;
                                        atom.orig_id  = mmdb2Atom[aNo]->serNum;
                                        atom.name     = mol_StrCopy( mmdb2Atom[aNo]->name, model ); // This assigns the allocated memory for the string to the model, while copying the string.
                                        mol_ResTypeConv ( mmdb2Residue[nRes]->name, &atom.res_type );
                                        atom.res_prop = mol_res_props[atom.res_type];
                                        
                                        if ( !( std::string ( mmdb2Atom[aNo]->altLoc ) == std::string ( "" ) ) &&
                                             !( std::string ( mmdb2Atom[aNo]->altLoc ) == std::string ( "A" ) ) )
                                            continue; // This is in keeping with the mol_DbPdbAtomProcLongChainId() function; it allows only the first alt-loc to be added.
                                        
                                        atom.res_seq  = mmdb2Residue[nRes]->seqNum;
                                        std::string ICodeHlp = std::string ( mmdb2Residue[nRes]->GetInsCode() );
                                        char ICode;
                                        if ( ICodeHlp == "" ) { ICode = ' '; } else { ICode = mmdb2Residue[nRes]->GetInsCode()[0]; }
                                        atom.insertion_code = ICode;
                                        atom.pos[0]   = mmdb2Atom[aNo]->x;
                                        atom.pos[1]   = mmdb2Atom[aNo]->y;
                                        atom.pos[2]   = mmdb2Atom[aNo]->z;
                                        atom.temp     = mmdb2Atom[aNo]->tempFactor;
                                        
                                        if (atom.name[0] != ' ') { mol_AtomTypeConv ( &atom.name[0], &atom.type ); }
                                        else                     { mol_AtomTypeConv ( &atom.name[1], &atom.type ); }
                                        
                                        if ( mmdb2Atom[aNo]->Het ) { atom.het = MOL_TRUE; }
                                        else                       { atom.het = MOL_FALSE; }
                                        atom.chain_id = mmdb2Chain[nCh]->GetChainID()[0];
                                        atom.long_chain_id = mol_StrCopy( " ", model ); // We do not have LongChainID, so empty
                                        
                                        mol_MolModelCurrStrucGet ( model, &struc );                 // Copies model->curr_struc to struc
                                        mol_StructureAtomAdd ( struc, mmdb2Atom[aNo]->Het, &atom ); // Saves the atom to the model (through struc)
                                    }
                                    else
                                    {
                                        std::cerr << "MMDB2 failed to read the atom " << aNo << " in residue " << nRes << " in chain " << nCh << " in file " << filename << ". The file is most likely corrupted. Terminating now..." << std::endl;
                                        exit                          ( -1 );
                                    }
                                }
                            }
                            else
                            {
                                std::cerr << "MMDB2 failed to read the residue " << nRes << " in chain " << nCh << " in file " << filename << ". The file is most likely corrupted. Terminating now..." << std::endl;
                                exit                          ( -1 );
                            }
                        }
                    }
                    else
                    {
                        std::cerr << "MMDB2 failed to read the chain " << nCh << " in file " << filename << ". The file is most likely corrupted. Terminating now..." << std::endl;
                        exit                          ( -1 );
                    }
                }
            }
            else
            {
                std::cerr << "MMDB2 failed to read the first model in file " << filename << ". The file is most likely corrupted. Terminating now..." << std::endl;
                exit                                  ( -1 );
            }
            
            //======================================== Build biomolecule from the filled in structs
            MolStructure *struc;
            mol_MolModelCurrStrucGet                  ( model, &struc ); // Fill in the struc structure
            mol_StructureChainsBuild                  ( struc, 1 );      // Build protein chains
            mol_StructureChainsBuild                  ( struc, 2 );      // Build solvent chains
            
            //======================================== Clean up
            delete mfile;
#else
            std::cerr << "ERROR: The suplied file has the mmCIF format, but Molmodel was not compiled with the Gemmi library, which is required for the mmCIF support. Please re-run the Molmodel cmake command with the -DUSE_GEMMI=TRUE -DGEMMI_PATH=/path/to/gemmi/include options (and then alse re-run the make install command). Terminating..." << std::endl;
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

