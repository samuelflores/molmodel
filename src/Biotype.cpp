/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Christopher Bruns                                                 *
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

// This is the one file where we want the handle instantiated
// #define DO_INSTANTIATE_BIOTYPE_PIMPL_HANDLE
// #include "molmodel/internal/Biotype.h"
// #undef DO_INSTANTIATE_BIOTYPE_PIMPL_HANDLE

#include "BiotypeRep.h"
// #include "SimTKcommon/internal/PrivateImplementation_Defs.h"

#include "molmodel/internal/common.h"
#include <map>
#include <algorithm>
#include <cctype>

using namespace std;

namespace SimTK {

// template class PIMPLImplementation<Biotype,BiotypeRep>; // explicit instantiation
// template class PIMPLHandle<Biotype,BiotypeRep>; // explicit instantiation

// Store biotype data in this compilation unit, rather than in the biotype class

static void initializePopularBiotypes();

static std::map<BiotypeIndex, Biotype> biotypesByIndex;
static std::map<TinkerBiotypeIndex, BiotypeIndex> biotypeIxsByTinkerIndex;
static BiotypeIndex nextUnusedBiotypeIndex(1);

class BiotypeKey {
public:
    BiotypeKey();
    BiotypeKey(const char* r, const char* a, Ordinality::Residue o)
        : residueName(regularizeString(r)), atomName(regularizeString(a)), ordinality(o)
    {}

    bool operator<(const BiotypeKey& other) const {
        if (residueName == other.residueName) {
            if (atomName == other.atomName)
                return (ordinality < other.ordinality);
            else
                return (atomName < other.atomName);
        }
        else {
            return (residueName < other.residueName);
        }
    }

    // convert to lower case and trim leading and trailing spaces
    static string regularizeString(const char* str)
    {
        size_t last = strlen(str);
        size_t first = 0;
        while (str[first] == ' ')
            first++;

        if (first == last)
            return "";

        last--;
        while (str[last] == ' ' && last > first)
            last--;

        string rs(last - first + 1, ' ');
        size_t idx = first;
        size_t jdx = 0;
        for (;idx <= last; idx++, jdx++) {
            char ch = str[idx];
            unsigned b = (unsigned(ch) - 65U) < 26U;
            rs[jdx] = ch + b * 32;
        }

        return rs;
    }

    String residueName;
    String atomName;
    Ordinality::Residue ordinality;
};

static bool tinkerIxExists(TinkerBiotypeIndex tIx) 
{
    return (biotypeIxsByTinkerIndex.find(tIx) != biotypeIxsByTinkerIndex.end());
}

static bool biotypeIxExists(BiotypeIndex bIx) 
{
    return (biotypesByIndex.find(bIx) != biotypesByIndex.end());
}

static std::map<BiotypeKey, BiotypeIndex> biotypeIxsByKey;

/* static */ BiotypeIndex Biotype::defineTinkerBiotype(
        TinkerBiotypeIndex tinkerBiotypeIndex, 
        const Element *element,
        int valence,
        const char* residueName, 
        const char* atomName, 
        Ordinality::Residue ordinality)
{
    assert(! Biotype::exists(residueName, atomName, ordinality) );

    BiotypeIndex biotypeIndex = nextUnusedBiotypeIndex;
    ++nextUnusedBiotypeIndex;
    assert(biotypesByIndex.find(biotypeIndex) == biotypesByIndex.end());

    BiotypeKey key(residueName, atomName, ordinality);

    biotypeIxsByKey[key] = biotypeIndex;
    biotypesByIndex[biotypeIndex] = Biotype(
        biotypeIndex,
        tinkerBiotypeIndex,
        element,
        valence,
        residueName,
        atomName,
        ordinality
        );

    if (tinkerBiotypeIndex.isValid()) {
        if (tinkerIxExists(tinkerBiotypeIndex))
		    throw std::range_error(
                String("TinkerBiotypeIndex already in use: ") + tinkerBiotypeIndex);
        assert(! tinkerIxExists(tinkerBiotypeIndex) );
        biotypeIxsByTinkerIndex[tinkerBiotypeIndex] = biotypeIndex;
        assert( tinkerIxExists(tinkerBiotypeIndex) );
    }

    assert( Biotype::exists(residueName, atomName, ordinality) );
    assert( biotypesByIndex.find(biotypeIndex) != biotypesByIndex.end() );

    return biotypeIndex;
}

static bool biotypeExists(const BiotypeKey& key) {
    return biotypeIxsByKey.find(key) != biotypeIxsByKey.end();
}

/* static */ bool Biotype::exists(const char* residueName, 
                   const char* atomName, 
                   Ordinality::Residue ordinality) {
    BiotypeKey key1(residueName, atomName, ordinality);
    BiotypeKey key2(residueName, atomName, Ordinality::Any);
    return (biotypeExists(key1) || biotypeExists(key2));
}

/* static */ bool Biotype::exists(BiotypeIndex biotypeIndex) {
    return biotypeIxExists(biotypeIndex);
}

////////////////
// Biotype::Impl //
////////////////

Biotype::Impl::Impl() {
}

Biotype::Impl::Impl(BiotypeIndex b,
                       TinkerBiotypeIndex tinkerBiotypeIndex, 
                       const Element* e,
                       int v,
                       const char* r, 
                       const char* a, 
                       Ordinality::Residue o)
     : biotypeIndex(b)
     , tinkerBiotypeIndexIfAny(tinkerBiotypeIndex)
     , element(e)
     , valence(v)
     , residueName(r)
     , atomName(a)
     , ordinality(o)
{}


///////////////
//  Biotype  //
///////////////

// These are all the constants of the Enumeration type Biotype.
///*static*/ const Biotype Biotype::Argon(Element::getBySymbol("Ar"), );
//
///*static*/ const Biotype Biotype::MethaneH;
///*static*/ const Biotype Biotype::MethaneC;
//
///*static*/ const Biotype Biotype::EthaneH;
///*static*/ const Biotype Biotype::EthaneC;
//
///*static*/ const Biotype Biotype::SerineN;
///*static*/ const Biotype Biotype::SerineHN;
///*static*/ const Biotype Biotype::SerineCA;
///*static*/ const Biotype Biotype::SerineHA;
///*static*/ const Biotype Biotype::SerineC;
///*static*/ const Biotype Biotype::SerineO;
///*static*/ const Biotype Biotype::SerineCB;
///*static*/ const Biotype Biotype::SerineHB;
///*static*/ const Biotype Biotype::SerineOG;
///*static*/ const Biotype Biotype::SerineHG;

Biotype::Biotype() 
    : impl(new Impl())
{}

Biotype::Biotype(BiotypeIndex b,
                 TinkerBiotypeIndex tinkerBiotypeIndex, 
                 const Element* e,
                 int v,
                 const char* r, 
                 const char* a, 
                 Ordinality::Residue o)
    : impl(new Impl(b, tinkerBiotypeIndex, e, v, r, a, o))
{}

const Element *  Biotype::getElement() const
{
    return impl->getElement();
}

int             Biotype::getValence() const
{
    return impl->getValence();
}

BiotypeIndex       Biotype::getIndex() const
{
    return impl->getIndex();
}

TinkerBiotypeIndex Biotype::getTinkerBiotypeIfAny() const
{
    return impl->getTinkerBiotypeIfAny();
}

Biotype& Biotype::setTinkerBiotypeIndex(TinkerBiotypeIndex tIx)
{
    impl->setTinkerBiotypeIndex(tIx);
    return *this;
}

const String& Biotype::getAtomName() const
{
    return impl->getAtomName();
}

const String& Biotype::getResidueName() const
{
    return impl->getResidueName();
}

Ordinality::Residue Biotype::getOrdinality() const
{
    return impl->getOrdinality();
}

std::ostream& Biotype::generateSelfCode(std::ostream& os) const
{
    return impl->generateSelfCode(os);
}

/* static */ const Biotype& Biotype::get(TinkerBiotypeIndex tinkerBiotypeIndex) {
    assert(tinkerIxExists(tinkerBiotypeIndex));
    BiotypeIndex id = biotypeIxsByTinkerIndex.find(tinkerBiotypeIndex)->second;
    assert(biotypeIxExists(id));
    return biotypesByIndex.find(id)->second;
}
 
/* static */ const Biotype& Biotype::get(const char* residueName, 
                          const char* atomName, 
                          Ordinality::Residue ordinality)
{
    assert(exists(residueName, atomName, ordinality));

    BiotypeKey key1(residueName, atomName, ordinality);
    BiotypeKey key2(residueName, atomName, Ordinality::Any);

    BiotypeIndex id;
    if (biotypeExists(key1))
        id = biotypeIxsByKey.find(key1)->second;
    else 
        id = biotypeIxsByKey.find(key2)->second;

    return get(id);
}

/* static */ const Biotype& Biotype::get(BiotypeIndex id)
{
    assert(biotypeIxExists(id));

    return biotypesByIndex.find(id)->second;
}

/* static */ Biotype& Biotype::upd(BiotypeIndex biotypeIndex) 
{
    assert(biotypeIxExists(biotypeIndex));
    return biotypesByIndex.find(biotypeIndex)->second;
}

/* static */ Biotype& Biotype::upd(TinkerBiotypeIndex tinkerBiotypeIndex)
{
    assert(tinkerIxExists(tinkerBiotypeIndex));
    BiotypeIndex id = biotypeIxsByTinkerIndex.find(tinkerBiotypeIndex)->second;
    return upd(id);
}

/* static */ Biotype& Biotype::upd(const char* residueName, 
                          const char* atomName, 
                          Ordinality::Residue ordinality)
{
    assert(exists(residueName, atomName, ordinality));

    BiotypeKey key1(residueName, atomName, ordinality);
    BiotypeKey key2(residueName, atomName, Ordinality::Any);

    BiotypeIndex id;
    if (biotypeExists(key1))
        id = biotypeIxsByKey.find(key1)->second;
    else 
        id = biotypeIxsByKey.find(key2)->second;

    assert(biotypeIxExists(id));

    return upd(id);
}


/* static */ std::ostream& Biotype::generateAllBiotypeCode(std::ostream& os) {
    std::map<BiotypeIndex, Biotype>::const_iterator biotypeI;
    for (biotypeI = biotypesByIndex.begin(); biotypeI != biotypesByIndex.end(); ++biotypeI) 
    {
        biotypeI->second.generateSelfCode(os);
    }

    return os;
}

 //Biotype::Biotype(const char* name,
//                 const Element& elt, int n, int chg, TinkerBiotypeIndex tid) 
//{
//    tinkerBiotypeIfAny = tid;
//    element            = elt;
//    valence            = n;
//    formalCharge       = chg;
//}

/* static */ const Biotype& Biotype::Argon() 
{
    initializePopularBiotypes();
    return Biotype::get("argon", "argon");
}

/* static */ const Biotype& Biotype::SerineO() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "O") )
        Biotype::defineBiotype(Element::getBySymbol("O"), 1, "serine", "O");
    return Biotype::get("serine", "O");
}

/* static */ const Biotype& Biotype::SerineC() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "C") )
        Biotype::defineBiotype(Element::getBySymbol("C"), 3, "serine", "C");
    return Biotype::get("serine", "C");
}

/* static */ const Biotype& Biotype::SerineHN() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "HN") )
        Biotype::defineBiotype(Element::getBySymbol("H"), 1, "serine", "HN");
    return Biotype::get("serine", "HN");
}

/* static */ const Biotype& Biotype::SerineHA() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "HA") )
        Biotype::defineBiotype(Element::getBySymbol("H"), 1, "serine", "HA");
    return Biotype::get("serine", "HA");
}

/* static */ const Biotype& Biotype::SerineHG() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "HG") )
        Biotype::defineBiotype(Element::getBySymbol("H"), 1, "serine", "HG");
    return Biotype::get("serine", "HG");
}

/* static */ const Biotype& Biotype::SerineOG() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "OG") )
        Biotype::defineBiotype(Element::getBySymbol("O"), 2, "serine", "OG");
    return Biotype::get("serine", "OG");
}

/* static */ const Biotype& Biotype::SerineN() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "N") )
        Biotype::defineBiotype(Element::getBySymbol("N"), 3, "serine", "N");
    return Biotype::get("serine", "N");
}

/* static */ const Biotype& Biotype::SerineCA() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("serine", "HA") )
        Biotype::defineBiotype(Element::getBySymbol("C"), 4, "serine", "CA");
    return Biotype::get("serine", "CA");
}

/* static */ const Biotype& Biotype::MethaneH() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("methane", "H") )
        Biotype::defineBiotype(Element::getBySymbol("H"), 1, "methane", "H");
    return Biotype::get("methane", "H");
}

/* static */ const Biotype& Biotype::MethaneC() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("methane", "C") )
        Biotype::defineBiotype(Element::getBySymbol("C"), 4, "methane", "C");
    return Biotype::get("methane", "C");
}

/* static */ const Biotype& Biotype::EthaneH() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("ethane", "H") )
        Biotype::defineBiotype(Element::getBySymbol("H"), 1, "ethane", "H");
    return Biotype::get("ethane", "H");
}

/* static */ const Biotype& Biotype::EthaneC() 
{
    initializePopularBiotypes();
    if (! Biotype::exists("ethane", "C") )
        Biotype::defineBiotype(Element::getBySymbol("C"), 4, "ethane", "C");
    return Biotype::get("ethane", "C");
}




/* static */ void Biotype::initializePopularBiotypes() 
{
    static bool popularBiotypesAreInitialized = false;
    if (popularBiotypesAreInitialized) return;

    if (! Biotype::exists("argon", "argon", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            InvalidTinkerBiotypeIndex
            , Element::getBySymbol("Ar")
            , 0
            , "argon"
            , "argon"
            , Ordinality::Any
            );

    if (! Biotype::exists("methane", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            InvalidTinkerBiotypeIndex
            , Element::getBySymbol("C")
            , 4
            , "methane"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("methane", "H", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            InvalidTinkerBiotypeIndex
            , Element::getBySymbol("H")
            , 1
            , "methane"
            , "H"
            , Ordinality::Any
            );

    if (! Biotype::exists("ethane", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            InvalidTinkerBiotypeIndex
            , Element::getBySymbol("C")
            , 4
            , "ethane"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("ethane", "H", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            InvalidTinkerBiotypeIndex
            , Element::getBySymbol("H")
            , 1
            , "ethane"
            , "H"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glycine", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1)
            , Element::getBySymbol("N")
            , 3
            , "Glycine"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glycine", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(2)
            , Element::getBySymbol("C")
            , 4
            , "Glycine"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glycine", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(3)
            , Element::getBySymbol("C")
            , 3
            , "Glycine"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glycine", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(4)
            , Element::getBySymbol("H")
            , 1
            , "Glycine"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glycine", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(5)
            , Element::getBySymbol("O")
            , 1
            , "Glycine"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glycine", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(6)
            , Element::getBySymbol("H")
            , 1
            , "Glycine"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Alanine", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(7)
            , Element::getBySymbol("N")
            , 3
            , "Alanine"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Alanine", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(8)
            , Element::getBySymbol("C")
            , 4
            , "Alanine"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Alanine", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(9)
            , Element::getBySymbol("C")
            , 3
            , "Alanine"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Alanine", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(10)
            , Element::getBySymbol("H")
            , 1
            , "Alanine"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Alanine", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(11)
            , Element::getBySymbol("O")
            , 1
            , "Alanine"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Alanine", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(12)
            , Element::getBySymbol("H")
            , 1
            , "Alanine"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Alanine", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(13)
            , Element::getBySymbol("C")
            , 4
            , "Alanine"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Alanine", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(14)
            , Element::getBySymbol("H")
            , 1
            , "Alanine"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Valine", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(15)
            , Element::getBySymbol("N")
            , 3
            , "Valine"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Valine", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(16)
            , Element::getBySymbol("C")
            , 4
            , "Valine"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Valine", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(17)
            , Element::getBySymbol("C")
            , 3
            , "Valine"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Valine", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(18)
            , Element::getBySymbol("H")
            , 1
            , "Valine"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Valine", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(19)
            , Element::getBySymbol("O")
            , 1
            , "Valine"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Valine", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(20)
            , Element::getBySymbol("H")
            , 1
            , "Valine"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Valine", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(21)
            , Element::getBySymbol("C")
            , 4
            , "Valine"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Valine", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(22)
            , Element::getBySymbol("H")
            , 1
            , "Valine"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Valine", "CG1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(23)
            , Element::getBySymbol("C")
            , 4
            , "Valine"
            , "CG1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Valine", "HG1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(24)
            , Element::getBySymbol("H")
            , 1
            , "Valine"
            , "HG1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Valine", "CG2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(25)
            , Element::getBySymbol("C")
            , 4
            , "Valine"
            , "CG2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Valine", "HG2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(26)
            , Element::getBySymbol("H")
            , 1
            , "Valine"
            , "HG2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Leucine", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(27)
            , Element::getBySymbol("N")
            , 3
            , "Leucine"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Leucine", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(28)
            , Element::getBySymbol("C")
            , 4
            , "Leucine"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Leucine", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(29)
            , Element::getBySymbol("C")
            , 3
            , "Leucine"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Leucine", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(30)
            , Element::getBySymbol("H")
            , 1
            , "Leucine"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Leucine", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(31)
            , Element::getBySymbol("O")
            , 1
            , "Leucine"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Leucine", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(32)
            , Element::getBySymbol("H")
            , 1
            , "Leucine"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Leucine", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(33)
            , Element::getBySymbol("C")
            , 4
            , "Leucine"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Leucine", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(34)
            , Element::getBySymbol("H")
            , 1
            , "Leucine"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Leucine", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(35)
            , Element::getBySymbol("C")
            , 4
            , "Leucine"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Leucine", "HG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(36)
            , Element::getBySymbol("H")
            , 1
            , "Leucine"
            , "HG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Leucine", "CD1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(37)
            , Element::getBySymbol("C")
            , 4
            , "Leucine"
            , "CD1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Leucine", "HD1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(38)
            , Element::getBySymbol("H")
            , 1
            , "Leucine"
            , "HD1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Leucine", "CD2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(39)
            , Element::getBySymbol("C")
            , 4
            , "Leucine"
            , "CD2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Leucine", "HD2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(40)
            , Element::getBySymbol("H")
            , 1
            , "Leucine"
            , "HD2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Isoleucine", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(41)
            , Element::getBySymbol("N")
            , 3
            , "Isoleucine"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Isoleucine", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(42)
            , Element::getBySymbol("C")
            , 4
            , "Isoleucine"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Isoleucine", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(43)
            , Element::getBySymbol("C")
            , 3
            , "Isoleucine"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Isoleucine", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(44)
            , Element::getBySymbol("H")
            , 1
            , "Isoleucine"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Isoleucine", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(45)
            , Element::getBySymbol("O")
            , 1
            , "Isoleucine"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Isoleucine", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(46)
            , Element::getBySymbol("H")
            , 1
            , "Isoleucine"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Isoleucine", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(47)
            , Element::getBySymbol("C")
            , 4
            , "Isoleucine"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Isoleucine", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(48)
            , Element::getBySymbol("H")
            , 1
            , "Isoleucine"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Isoleucine", "CG1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(49)
            , Element::getBySymbol("C")
            , 4
            , "Isoleucine"
            , "CG1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Isoleucine", "HG1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(50)
            , Element::getBySymbol("H")
            , 1
            , "Isoleucine"
            , "HG1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Isoleucine", "CG2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(51)
            , Element::getBySymbol("C")
            , 4
            , "Isoleucine"
            , "CG2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Isoleucine", "HG2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(52)
            , Element::getBySymbol("H")
            , 1
            , "Isoleucine"
            , "HG2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Isoleucine", "CD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(53)
            , Element::getBySymbol("C")
            , 4
            , "Isoleucine"
            , "CD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Isoleucine", "HD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(54)
            , Element::getBySymbol("H")
            , 1
            , "Isoleucine"
            , "HD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Serine", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(55)
            , Element::getBySymbol("N")
            , 3
            , "Serine"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Serine", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(56)
            , Element::getBySymbol("C")
            , 4
            , "Serine"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Serine", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(57)
            , Element::getBySymbol("C")
            , 3
            , "Serine"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Serine", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(58)
            , Element::getBySymbol("H")
            , 1
            , "Serine"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Serine", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(59)
            , Element::getBySymbol("O")
            , 1
            , "Serine"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Serine", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(60)
            , Element::getBySymbol("H")
            , 1
            , "Serine"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Serine", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(61)
            , Element::getBySymbol("C")
            , 4
            , "Serine"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Serine", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(62)
            , Element::getBySymbol("H")
            , 1
            , "Serine"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Serine", "OG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(63)
            , Element::getBySymbol("O")
            , 2
            , "Serine"
            , "OG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Serine", "HG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(64)
            , Element::getBySymbol("H")
            , 1
            , "Serine"
            , "HG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Threonine", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(65)
            , Element::getBySymbol("N")
            , 3
            , "Threonine"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Threonine", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(66)
            , Element::getBySymbol("C")
            , 4
            , "Threonine"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Threonine", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(67)
            , Element::getBySymbol("C")
            , 3
            , "Threonine"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Threonine", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(68)
            , Element::getBySymbol("H")
            , 1
            , "Threonine"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Threonine", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(69)
            , Element::getBySymbol("O")
            , 1
            , "Threonine"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Threonine", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(70)
            , Element::getBySymbol("H")
            , 1
            , "Threonine"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Threonine", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(71)
            , Element::getBySymbol("C")
            , 4
            , "Threonine"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Threonine", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(72)
            , Element::getBySymbol("H")
            , 1
            , "Threonine"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Threonine", "OG1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(73)
            , Element::getBySymbol("O")
            , 2
            , "Threonine"
            , "OG1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Threonine", "HG1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(74)
            , Element::getBySymbol("H")
            , 1
            , "Threonine"
            , "HG1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Threonine", "CG2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(75)
            , Element::getBySymbol("C")
            , 4
            , "Threonine"
            , "CG2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Threonine", "HG2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(76)
            , Element::getBySymbol("H")
            , 1
            , "Threonine"
            , "HG2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cysteine (-SH)", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(77)
            , Element::getBySymbol("N")
            , 3
            , "Cysteine (-SH)"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cysteine (-SH)", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(78)
            , Element::getBySymbol("C")
            , 4
            , "Cysteine (-SH)"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cysteine (-SH)", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(79)
            , Element::getBySymbol("C")
            , 3
            , "Cysteine (-SH)"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cysteine (-SH)", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(80)
            , Element::getBySymbol("H")
            , 1
            , "Cysteine (-SH)"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cysteine (-SH)", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(81)
            , Element::getBySymbol("O")
            , 1
            , "Cysteine (-SH)"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cysteine (-SH)", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(82)
            , Element::getBySymbol("H")
            , 1
            , "Cysteine (-SH)"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cysteine (-SH)", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(83)
            , Element::getBySymbol("C")
            , 4
            , "Cysteine (-SH)"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cysteine (-SH)", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(84)
            , Element::getBySymbol("H")
            , 1
            , "Cysteine (-SH)"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cysteine (-SH)", "SG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(85)
            , Element::getBySymbol("S")
            , 2
            , "Cysteine (-SH)"
            , "SG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cysteine (-SH)", "HG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(86)
            , Element::getBySymbol("H")
            , 1
            , "Cysteine (-SH)"
            , "HG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cystine (-SS-)", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(87)
            , Element::getBySymbol("N")
            , 3
            , "Cystine (-SS-)"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cystine (-SS-)", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(88)
            , Element::getBySymbol("C")
            , 4
            , "Cystine (-SS-)"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cystine (-SS-)", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(89)
            , Element::getBySymbol("C")
            , 3
            , "Cystine (-SS-)"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cystine (-SS-)", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(90)
            , Element::getBySymbol("H")
            , 1
            , "Cystine (-SS-)"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cystine (-SS-)", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(91)
            , Element::getBySymbol("O")
            , 1
            , "Cystine (-SS-)"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cystine (-SS-)", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(92)
            , Element::getBySymbol("H")
            , 1
            , "Cystine (-SS-)"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cystine (-SS-)", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(93)
            , Element::getBySymbol("C")
            , 4
            , "Cystine (-SS-)"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cystine (-SS-)", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(94)
            , Element::getBySymbol("H")
            , 1
            , "Cystine (-SS-)"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cystine (-SS-)", "SG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(95)
            , Element::getBySymbol("S")
            , 2
            , "Cystine (-SS-)"
            , "SG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Proline", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(96)
            , Element::getBySymbol("N")
            , 3
            , "Proline"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Proline", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(97)
            , Element::getBySymbol("C")
            , 4
            , "Proline"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Proline", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(98)
            , Element::getBySymbol("C")
            , 3
            , "Proline"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Proline", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(99)
            , Element::getBySymbol("O")
            , 1
            , "Proline"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Proline", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(100)
            , Element::getBySymbol("H")
            , 1
            , "Proline"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Proline", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(101)
            , Element::getBySymbol("C")
            , 4
            , "Proline"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Proline", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(102)
            , Element::getBySymbol("H")
            , 1
            , "Proline"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Proline", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(103)
            , Element::getBySymbol("C")
            , 4
            , "Proline"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Proline", "HG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(104)
            , Element::getBySymbol("H")
            , 1
            , "Proline"
            , "HG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Proline", "CD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(105)
            , Element::getBySymbol("C")
            , 4
            , "Proline"
            , "CD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Proline", "HD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(106)
            , Element::getBySymbol("H")
            , 1
            , "Proline"
            , "HD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phenylalanine", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(107)
            , Element::getBySymbol("N")
            , 3
            , "Phenylalanine"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phenylalanine", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(108)
            , Element::getBySymbol("C")
            , 4
            , "Phenylalanine"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phenylalanine", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(109)
            , Element::getBySymbol("C")
            , 3
            , "Phenylalanine"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phenylalanine", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(110)
            , Element::getBySymbol("H")
            , 1
            , "Phenylalanine"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phenylalanine", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(111)
            , Element::getBySymbol("O")
            , 1
            , "Phenylalanine"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phenylalanine", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(112)
            , Element::getBySymbol("H")
            , 1
            , "Phenylalanine"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phenylalanine", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(113)
            , Element::getBySymbol("C")
            , 4
            , "Phenylalanine"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phenylalanine", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(114)
            , Element::getBySymbol("H")
            , 1
            , "Phenylalanine"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phenylalanine", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(115)
            , Element::getBySymbol("C")
            , 3
            , "Phenylalanine"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phenylalanine", "CD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(116)
            , Element::getBySymbol("C")
            , 3
            , "Phenylalanine"
            , "CD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phenylalanine", "HD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(117)
            , Element::getBySymbol("H")
            , 1
            , "Phenylalanine"
            , "HD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phenylalanine", "CE", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(118)
            , Element::getBySymbol("C")
            , 3
            , "Phenylalanine"
            , "CE"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phenylalanine", "HE", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(119)
            , Element::getBySymbol("H")
            , 1
            , "Phenylalanine"
            , "HE"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phenylalanine", "CZ", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(120)
            , Element::getBySymbol("C")
            , 3
            , "Phenylalanine"
            , "CZ"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phenylalanine", "HZ", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(121)
            , Element::getBySymbol("H")
            , 1
            , "Phenylalanine"
            , "HZ"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(122)
            , Element::getBySymbol("N")
            , 3
            , "Tyrosine"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(123)
            , Element::getBySymbol("C")
            , 4
            , "Tyrosine"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(124)
            , Element::getBySymbol("C")
            , 3
            , "Tyrosine"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(125)
            , Element::getBySymbol("H")
            , 1
            , "Tyrosine"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(126)
            , Element::getBySymbol("O")
            , 1
            , "Tyrosine"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(127)
            , Element::getBySymbol("H")
            , 1
            , "Tyrosine"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(128)
            , Element::getBySymbol("C")
            , 4
            , "Tyrosine"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(129)
            , Element::getBySymbol("H")
            , 1
            , "Tyrosine"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(130)
            , Element::getBySymbol("C")
            , 3
            , "Tyrosine"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "CD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(131)
            , Element::getBySymbol("C")
            , 3
            , "Tyrosine"
            , "CD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "HD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(132)
            , Element::getBySymbol("H")
            , 1
            , "Tyrosine"
            , "HD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "CE", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(133)
            , Element::getBySymbol("C")
            , 3
            , "Tyrosine"
            , "CE"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "HE", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(134)
            , Element::getBySymbol("H")
            , 1
            , "Tyrosine"
            , "HE"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "CZ", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(135)
            , Element::getBySymbol("C")
            , 3
            , "Tyrosine"
            , "CZ"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "OH", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(136)
            , Element::getBySymbol("O")
            , 2
            , "Tyrosine"
            , "OH"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tyrosine", "HH", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(137)
            , Element::getBySymbol("H")
            , 1
            , "Tyrosine"
            , "HH"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(138)
            , Element::getBySymbol("N")
            , 3
            , "Tryptophan"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(139)
            , Element::getBySymbol("C")
            , 4
            , "Tryptophan"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(140)
            , Element::getBySymbol("C")
            , 3
            , "Tryptophan"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(141)
            , Element::getBySymbol("H")
            , 1
            , "Tryptophan"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(142)
            , Element::getBySymbol("O")
            , 1
            , "Tryptophan"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(143)
            , Element::getBySymbol("H")
            , 1
            , "Tryptophan"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(144)
            , Element::getBySymbol("C")
            , 4
            , "Tryptophan"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(145)
            , Element::getBySymbol("H")
            , 1
            , "Tryptophan"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(146)
            , Element::getBySymbol("C")
            , 3
            , "Tryptophan"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "CD1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(147)
            , Element::getBySymbol("C")
            , 3
            , "Tryptophan"
            , "CD1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "HD1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(148)
            , Element::getBySymbol("H")
            , 1
            , "Tryptophan"
            , "HD1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "CD2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(149)
            , Element::getBySymbol("C")
            , 3
            , "Tryptophan"
            , "CD2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "NE1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(150)
            , Element::getBySymbol("N")
            , 3
            , "Tryptophan"
            , "NE1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "HE1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(151)
            , Element::getBySymbol("H")
            , 1
            , "Tryptophan"
            , "HE1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "CE2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(152)
            , Element::getBySymbol("C")
            , 3
            , "Tryptophan"
            , "CE2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "CE3", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(153)
            , Element::getBySymbol("C")
            , 3
            , "Tryptophan"
            , "CE3"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "HE3", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(154)
            , Element::getBySymbol("H")
            , 1
            , "Tryptophan"
            , "HE3"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "CZ2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(155)
            , Element::getBySymbol("C")
            , 3
            , "Tryptophan"
            , "CZ2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "HZ2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(156)
            , Element::getBySymbol("H")
            , 1
            , "Tryptophan"
            , "HZ2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "CZ3", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(157)
            , Element::getBySymbol("C")
            , 3
            , "Tryptophan"
            , "CZ3"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "HZ3", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(158)
            , Element::getBySymbol("H")
            , 1
            , "Tryptophan"
            , "HZ3"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "CH2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(159)
            , Element::getBySymbol("C")
            , 3
            , "Tryptophan"
            , "CH2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Tryptophan", "HH2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(160)
            , Element::getBySymbol("H")
            , 1
            , "Tryptophan"
            , "HH2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(161)
            , Element::getBySymbol("N")
            , 3
            , "Histidine (+)"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(162)
            , Element::getBySymbol("C")
            , 4
            , "Histidine (+)"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(163)
            , Element::getBySymbol("C")
            , 3
            , "Histidine (+)"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(164)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (+)"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(165)
            , Element::getBySymbol("O")
            , 1
            , "Histidine (+)"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(166)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (+)"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(167)
            , Element::getBySymbol("C")
            , 4
            , "Histidine (+)"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(168)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (+)"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(169)
            , Element::getBySymbol("C")
            , 3
            , "Histidine (+)"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "ND1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(170)
            , Element::getBySymbol("N")
            , 3
            , "Histidine (+)"
            , "ND1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "HD1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(171)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (+)"
            , "HD1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "CD2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(172)
            , Element::getBySymbol("C")
            , 3
            , "Histidine (+)"
            , "CD2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "HD2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(173)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (+)"
            , "HD2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "CE1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(174)
            , Element::getBySymbol("C")
            , 3
            , "Histidine (+)"
            , "CE1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "HE1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(175)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (+)"
            , "HE1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "NE2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(176)
            , Element::getBySymbol("N")
            , 3
            , "Histidine (+)"
            , "NE2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (+)", "HE2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(177)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (+)"
            , "HE2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(178)
            , Element::getBySymbol("N")
            , 3
            , "Histidine (HD)"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(179)
            , Element::getBySymbol("C")
            , 4
            , "Histidine (HD)"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(180)
            , Element::getBySymbol("C")
            , 3
            , "Histidine (HD)"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(181)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (HD)"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(182)
            , Element::getBySymbol("O")
            , 1
            , "Histidine (HD)"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(183)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (HD)"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(184)
            , Element::getBySymbol("C")
            , 4
            , "Histidine (HD)"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(185)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (HD)"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(186)
            , Element::getBySymbol("C")
            , 3
            , "Histidine (HD)"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "ND1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(187)
            , Element::getBySymbol("N")
            , 3
            , "Histidine (HD)"
            , "ND1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "HD1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(188)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (HD)"
            , "HD1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "CD2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(189)
            , Element::getBySymbol("C")
            , 3
            , "Histidine (HD)"
            , "CD2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "HD2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(190)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (HD)"
            , "HD2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "CE1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(191)
            , Element::getBySymbol("C")
            , 3
            , "Histidine (HD)"
            , "CE1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "HE1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(192)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (HD)"
            , "HE1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HD)", "NE2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(193)
            , Element::getBySymbol("N")
            , 2
            , "Histidine (HD)"
            , "NE2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(194)
            , Element::getBySymbol("N")
            , 3
            , "Histidine (HE)"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(195)
            , Element::getBySymbol("C")
            , 4
            , "Histidine (HE)"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(196)
            , Element::getBySymbol("C")
            , 3
            , "Histidine (HE)"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(197)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (HE)"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(198)
            , Element::getBySymbol("O")
            , 1
            , "Histidine (HE)"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(199)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (HE)"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(200)
            , Element::getBySymbol("C")
            , 4
            , "Histidine (HE)"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(201)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (HE)"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(202)
            , Element::getBySymbol("C")
            , 3
            , "Histidine (HE)"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "ND1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(203)
            , Element::getBySymbol("N")
            , 2
            , "Histidine (HE)"
            , "ND1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "CD2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(204)
            , Element::getBySymbol("C")
            , 3
            , "Histidine (HE)"
            , "CD2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "HD2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(205)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (HE)"
            , "HD2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "CE1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(206)
            , Element::getBySymbol("C")
            , 3
            , "Histidine (HE)"
            , "CE1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "HE1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(207)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (HE)"
            , "HE1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "NE2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(208)
            , Element::getBySymbol("N")
            , 3
            , "Histidine (HE)"
            , "NE2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Histidine (HE)", "HE2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(209)
            , Element::getBySymbol("H")
            , 1
            , "Histidine (HE)"
            , "HE2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Aspartic Acid", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(210)
            , Element::getBySymbol("N")
            , 3
            , "Aspartic Acid"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Aspartic Acid", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(211)
            , Element::getBySymbol("C")
            , 4
            , "Aspartic Acid"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Aspartic Acid", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(212)
            , Element::getBySymbol("C")
            , 3
            , "Aspartic Acid"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Aspartic Acid", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(213)
            , Element::getBySymbol("H")
            , 1
            , "Aspartic Acid"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Aspartic Acid", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(214)
            , Element::getBySymbol("O")
            , 1
            , "Aspartic Acid"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Aspartic Acid", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(215)
            , Element::getBySymbol("H")
            , 1
            , "Aspartic Acid"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Aspartic Acid", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(216)
            , Element::getBySymbol("C")
            , 4
            , "Aspartic Acid"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Aspartic Acid", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(217)
            , Element::getBySymbol("H")
            , 1
            , "Aspartic Acid"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Aspartic Acid", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(218)
            , Element::getBySymbol("C")
            , 3
            , "Aspartic Acid"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Aspartic Acid", "OD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(219)
            , Element::getBySymbol("O")
            , 1
            , "Aspartic Acid"
            , "OD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Asparagine", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(220)
            , Element::getBySymbol("N")
            , 3
            , "Asparagine"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Asparagine", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(221)
            , Element::getBySymbol("C")
            , 4
            , "Asparagine"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Asparagine", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(222)
            , Element::getBySymbol("C")
            , 3
            , "Asparagine"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Asparagine", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(223)
            , Element::getBySymbol("H")
            , 1
            , "Asparagine"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Asparagine", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(224)
            , Element::getBySymbol("O")
            , 1
            , "Asparagine"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Asparagine", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(225)
            , Element::getBySymbol("H")
            , 1
            , "Asparagine"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Asparagine", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(226)
            , Element::getBySymbol("C")
            , 4
            , "Asparagine"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Asparagine", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(227)
            , Element::getBySymbol("H")
            , 1
            , "Asparagine"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Asparagine", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(228)
            , Element::getBySymbol("C")
            , 3
            , "Asparagine"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Asparagine", "OD1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(229)
            , Element::getBySymbol("O")
            , 1
            , "Asparagine"
            , "OD1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Asparagine", "ND2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(230)
            , Element::getBySymbol("N")
            , 3
            , "Asparagine"
            , "ND2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Asparagine", "HD2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(231)
            , Element::getBySymbol("H")
            , 1
            , "Asparagine"
            , "HD2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamic Acid", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(232)
            , Element::getBySymbol("N")
            , 3
            , "Glutamic Acid"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamic Acid", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(233)
            , Element::getBySymbol("C")
            , 4
            , "Glutamic Acid"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamic Acid", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(234)
            , Element::getBySymbol("C")
            , 3
            , "Glutamic Acid"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamic Acid", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(235)
            , Element::getBySymbol("H")
            , 1
            , "Glutamic Acid"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamic Acid", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(236)
            , Element::getBySymbol("O")
            , 1
            , "Glutamic Acid"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamic Acid", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(237)
            , Element::getBySymbol("H")
            , 1
            , "Glutamic Acid"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamic Acid", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(238)
            , Element::getBySymbol("C")
            , 4
            , "Glutamic Acid"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamic Acid", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(239)
            , Element::getBySymbol("H")
            , 1
            , "Glutamic Acid"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamic Acid", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(240)
            , Element::getBySymbol("C")
            , 4
            , "Glutamic Acid"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamic Acid", "HG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(241)
            , Element::getBySymbol("H")
            , 1
            , "Glutamic Acid"
            , "HG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamic Acid", "CD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(242)
            , Element::getBySymbol("C")
            , 3
            , "Glutamic Acid"
            , "CD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamic Acid", "OE", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(243)
            , Element::getBySymbol("O")
            , 1
            , "Glutamic Acid"
            , "OE"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamine", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(244)
            , Element::getBySymbol("N")
            , 3
            , "Glutamine"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamine", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(245)
            , Element::getBySymbol("C")
            , 4
            , "Glutamine"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamine", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(246)
            , Element::getBySymbol("C")
            , 3
            , "Glutamine"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamine", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(247)
            , Element::getBySymbol("H")
            , 1
            , "Glutamine"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamine", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(248)
            , Element::getBySymbol("O")
            , 1
            , "Glutamine"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamine", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(249)
            , Element::getBySymbol("H")
            , 1
            , "Glutamine"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamine", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(250)
            , Element::getBySymbol("C")
            , 4
            , "Glutamine"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamine", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(251)
            , Element::getBySymbol("H")
            , 1
            , "Glutamine"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamine", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(252)
            , Element::getBySymbol("C")
            , 4
            , "Glutamine"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamine", "HG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(253)
            , Element::getBySymbol("H")
            , 1
            , "Glutamine"
            , "HG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamine", "CD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(254)
            , Element::getBySymbol("C")
            , 3
            , "Glutamine"
            , "CD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamine", "OE1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(255)
            , Element::getBySymbol("O")
            , 1
            , "Glutamine"
            , "OE1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamine", "NE2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(256)
            , Element::getBySymbol("N")
            , 3
            , "Glutamine"
            , "NE2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Glutamine", "HE2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(257)
            , Element::getBySymbol("H")
            , 1
            , "Glutamine"
            , "HE2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Methionine", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(258)
            , Element::getBySymbol("N")
            , 3
            , "Methionine"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Methionine", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(259)
            , Element::getBySymbol("C")
            , 4
            , "Methionine"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Methionine", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(260)
            , Element::getBySymbol("C")
            , 3
            , "Methionine"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Methionine", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(261)
            , Element::getBySymbol("H")
            , 1
            , "Methionine"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Methionine", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(262)
            , Element::getBySymbol("O")
            , 1
            , "Methionine"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Methionine", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(263)
            , Element::getBySymbol("H")
            , 1
            , "Methionine"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Methionine", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(264)
            , Element::getBySymbol("C")
            , 4
            , "Methionine"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Methionine", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(265)
            , Element::getBySymbol("H")
            , 1
            , "Methionine"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Methionine", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(266)
            , Element::getBySymbol("C")
            , 4
            , "Methionine"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Methionine", "HG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(267)
            , Element::getBySymbol("H")
            , 1
            , "Methionine"
            , "HG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Methionine", "SD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(268)
            , Element::getBySymbol("S")
            , 2
            , "Methionine"
            , "SD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Methionine", "CE", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(269)
            , Element::getBySymbol("C")
            , 4
            , "Methionine"
            , "CE"
            , Ordinality::Any
            );

    if (! Biotype::exists("Methionine", "HE", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(270)
            , Element::getBySymbol("H")
            , 1
            , "Methionine"
            , "HE"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(271)
            , Element::getBySymbol("N")
            , 3
            , "Lysine"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(272)
            , Element::getBySymbol("C")
            , 4
            , "Lysine"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(273)
            , Element::getBySymbol("C")
            , 3
            , "Lysine"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(274)
            , Element::getBySymbol("H")
            , 1
            , "Lysine"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(275)
            , Element::getBySymbol("O")
            , 1
            , "Lysine"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(276)
            , Element::getBySymbol("H")
            , 1
            , "Lysine"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(277)
            , Element::getBySymbol("C")
            , 4
            , "Lysine"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(278)
            , Element::getBySymbol("H")
            , 1
            , "Lysine"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(279)
            , Element::getBySymbol("C")
            , 4
            , "Lysine"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "HG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(280)
            , Element::getBySymbol("H")
            , 1
            , "Lysine"
            , "HG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "CD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(281)
            , Element::getBySymbol("C")
            , 4
            , "Lysine"
            , "CD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "HD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(282)
            , Element::getBySymbol("H")
            , 1
            , "Lysine"
            , "HD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "CE", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(283)
            , Element::getBySymbol("C")
            , 4
            , "Lysine"
            , "CE"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "HE", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(284)
            , Element::getBySymbol("H")
            , 1
            , "Lysine"
            , "HE"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "NZ", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(285)
            , Element::getBySymbol("N")
            , 4
            , "Lysine"
            , "NZ"
            , Ordinality::Any
            );

    if (! Biotype::exists("Lysine", "HZ", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(286)
            , Element::getBySymbol("H")
            , 1
            , "Lysine"
            , "HZ"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(287)
            , Element::getBySymbol("N")
            , 3
            , "Arginine"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(288)
            , Element::getBySymbol("C")
            , 4
            , "Arginine"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(289)
            , Element::getBySymbol("C")
            , 3
            , "Arginine"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(290)
            , Element::getBySymbol("H")
            , 1
            , "Arginine"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(291)
            , Element::getBySymbol("O")
            , 1
            , "Arginine"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(292)
            , Element::getBySymbol("H")
            , 1
            , "Arginine"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(293)
            , Element::getBySymbol("C")
            , 4
            , "Arginine"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(294)
            , Element::getBySymbol("H")
            , 1
            , "Arginine"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(295)
            , Element::getBySymbol("C")
            , 4
            , "Arginine"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "HG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(296)
            , Element::getBySymbol("H")
            , 1
            , "Arginine"
            , "HG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "CD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(297)
            , Element::getBySymbol("C")
            , 4
            , "Arginine"
            , "CD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "HD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(298)
            , Element::getBySymbol("H")
            , 1
            , "Arginine"
            , "HD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "NE", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(299)
            , Element::getBySymbol("N")
            , 3
            , "Arginine"
            , "NE"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "HE", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(300)
            , Element::getBySymbol("H")
            , 1
            , "Arginine"
            , "HE"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "CZ", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(301)
            , Element::getBySymbol("C")
            , 3
            , "Arginine"
            , "CZ"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "NH", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(302)
            , Element::getBySymbol("N")
            , 3
            , "Arginine"
            , "NH"
            , Ordinality::Any
            );

    if (! Biotype::exists("Arginine", "HH", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(303)
            , Element::getBySymbol("H")
            , 1
            , "Arginine"
            , "HH"
            , Ordinality::Any
            );

    if (! Biotype::exists("Ornithine", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(304)
            , Element::getBySymbol("N")
            , 3
            , "Ornithine"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Ornithine", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(305)
            , Element::getBySymbol("C")
            , 4
            , "Ornithine"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Ornithine", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(306)
            , Element::getBySymbol("C")
            , 3
            , "Ornithine"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Ornithine", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(307)
            , Element::getBySymbol("H")
            , 1
            , "Ornithine"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Ornithine", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(308)
            , Element::getBySymbol("O")
            , 1
            , "Ornithine"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Ornithine", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(309)
            , Element::getBySymbol("H")
            , 1
            , "Ornithine"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Ornithine", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(310)
            , Element::getBySymbol("C")
            , 4
            , "Ornithine"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Ornithine", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(311)
            , Element::getBySymbol("H")
            , 1
            , "Ornithine"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Ornithine", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(312)
            , Element::getBySymbol("C")
            , 4
            , "Ornithine"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Ornithine", "HG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(313)
            , Element::getBySymbol("H")
            , 1
            , "Ornithine"
            , "HG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Ornithine", "CD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(314)
            , Element::getBySymbol("C")
            , 4
            , "Ornithine"
            , "CD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Ornithine", "HD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(315)
            , Element::getBySymbol("H")
            , 1
            , "Ornithine"
            , "HD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Ornithine", "NE", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(316)
            , Element::getBySymbol("N")
            , 4
            , "Ornithine"
            , "NE"
            , Ordinality::Any
            );

    if (! Biotype::exists("Ornithine", "HE", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(317)
            , Element::getBySymbol("H")
            , 1
            , "Ornithine"
            , "HE"
            , Ordinality::Any
            );

    if (! Biotype::exists("MethylAlanine (AIB)", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(318)
            , Element::getBySymbol("N")
            , 3
            , "MethylAlanine (AIB)"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("MethylAlanine (AIB)", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(319)
            , Element::getBySymbol("C")
            , 4
            , "MethylAlanine (AIB)"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("MethylAlanine (AIB)", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(320)
            , Element::getBySymbol("C")
            , 3
            , "MethylAlanine (AIB)"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("MethylAlanine (AIB)", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(321)
            , Element::getBySymbol("H")
            , 1
            , "MethylAlanine (AIB)"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("MethylAlanine (AIB)", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(322)
            , Element::getBySymbol("O")
            , 1
            , "MethylAlanine (AIB)"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("MethylAlanine (AIB)", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(323)
            , Element::getBySymbol("C")
            , 4
            , "MethylAlanine (AIB)"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("MethylAlanine (AIB)", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(324)
            , Element::getBySymbol("H")
            , 1
            , "MethylAlanine (AIB)"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Pyroglutamic Acid", "N", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(325)
            , Element::getBySymbol("N")
            , 3
            , "Pyroglutamic Acid"
            , "N"
            , Ordinality::Any
            );

    if (! Biotype::exists("Pyroglutamic Acid", "CA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(326)
            , Element::getBySymbol("C")
            , 4
            , "Pyroglutamic Acid"
            , "CA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Pyroglutamic Acid", "C", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(327)
            , Element::getBySymbol("C")
            , 3
            , "Pyroglutamic Acid"
            , "C"
            , Ordinality::Any
            );

    if (! Biotype::exists("Pyroglutamic Acid", "HN", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(328)
            , Element::getBySymbol("H")
            , 1
            , "Pyroglutamic Acid"
            , "HN"
            , Ordinality::Any
            );

    if (! Biotype::exists("Pyroglutamic Acid", "O", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(329)
            , Element::getBySymbol("O")
            , 1
            , "Pyroglutamic Acid"
            , "O"
            , Ordinality::Any
            );

    if (! Biotype::exists("Pyroglutamic Acid", "HA", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(330)
            , Element::getBySymbol("H")
            , 1
            , "Pyroglutamic Acid"
            , "HA"
            , Ordinality::Any
            );

    if (! Biotype::exists("Pyroglutamic Acid", "CB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(331)
            , Element::getBySymbol("C")
            , 4
            , "Pyroglutamic Acid"
            , "CB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Pyroglutamic Acid", "HB", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(332)
            , Element::getBySymbol("H")
            , 1
            , "Pyroglutamic Acid"
            , "HB"
            , Ordinality::Any
            );

    if (! Biotype::exists("Pyroglutamic Acid", "CG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(333)
            , Element::getBySymbol("C")
            , 4
            , "Pyroglutamic Acid"
            , "CG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Pyroglutamic Acid", "HG", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(334)
            , Element::getBySymbol("H")
            , 1
            , "Pyroglutamic Acid"
            , "HG"
            , Ordinality::Any
            );

    if (! Biotype::exists("Pyroglutamic Acid", "CD", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(335)
            , Element::getBySymbol("C")
            , 3
            , "Pyroglutamic Acid"
            , "CD"
            , Ordinality::Any
            );

    if (! Biotype::exists("Pyroglutamic Acid", "OE", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(336)
            , Element::getBySymbol("O")
            , 1
            , "Pyroglutamic Acid"
            , "OE"
            , Ordinality::Any
            );

    if (! Biotype::exists("Acetyl", "CH3", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(340)
            , Element::getBySymbol("C")
            , 4
            , "Acetyl"
            , "CH3"
            , Ordinality::Initial
            );

    if (! Biotype::exists("Acetyl", "H", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(341)
            , Element::getBySymbol("H")
            , 1
            , "Acetyl"
            , "H"
            , Ordinality::Initial
            );

    if (! Biotype::exists("Acetyl", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(342)
            , Element::getBySymbol("C")
            , 3
            , "Acetyl"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("Acetyl", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(343)
            , Element::getBySymbol("O")
            , 1
            , "Acetyl"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("Amide", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(344)
            , Element::getBySymbol("N")
            , 3
            , "Amide"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("Amide", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(345)
            , Element::getBySymbol("H")
            , 1
            , "Amide"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("N-MeAmide", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(346)
            , Element::getBySymbol("N")
            , 3
            , "N-MeAmide"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("N-MeAmide", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(347)
            , Element::getBySymbol("H")
            , 1
            , "N-MeAmide"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("N-MeAmide", "CH3", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(348)
            , Element::getBySymbol("C")
            , 4
            , "N-MeAmide"
            , "CH3"
            , Ordinality::Final
            );

    if (! Biotype::exists("N-MeAmide", "H", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(349)
            , Element::getBySymbol("H")
            , 1
            , "N-MeAmide"
            , "H"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLY", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(350)
            , Element::getBySymbol("N")
            , 4
            , "GLY"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLY", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(351)
            , Element::getBySymbol("C")
            , 4
            , "GLY"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLY", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(352)
            , Element::getBySymbol("C")
            , 3
            , "GLY"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLY", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(353)
            , Element::getBySymbol("H")
            , 1
            , "GLY"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLY", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(354)
            , Element::getBySymbol("O")
            , 1
            , "GLY"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLY", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(355)
            , Element::getBySymbol("H")
            , 1
            , "GLY"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ALA", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(356)
            , Element::getBySymbol("N")
            , 4
            , "ALA"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ALA", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(357)
            , Element::getBySymbol("C")
            , 4
            , "ALA"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ALA", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(358)
            , Element::getBySymbol("C")
            , 3
            , "ALA"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ALA", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(359)
            , Element::getBySymbol("H")
            , 1
            , "ALA"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ALA", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(360)
            , Element::getBySymbol("O")
            , 1
            , "ALA"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ALA", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(361)
            , Element::getBySymbol("H")
            , 1
            , "ALA"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("VAL", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(362)
            , Element::getBySymbol("N")
            , 4
            , "VAL"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("VAL", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(363)
            , Element::getBySymbol("C")
            , 4
            , "VAL"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("VAL", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(364)
            , Element::getBySymbol("C")
            , 3
            , "VAL"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("VAL", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(365)
            , Element::getBySymbol("H")
            , 1
            , "VAL"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("VAL", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(366)
            , Element::getBySymbol("O")
            , 1
            , "VAL"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("VAL", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(367)
            , Element::getBySymbol("H")
            , 1
            , "VAL"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("LEU", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(368)
            , Element::getBySymbol("N")
            , 4
            , "LEU"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("LEU", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(369)
            , Element::getBySymbol("C")
            , 4
            , "LEU"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("LEU", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(370)
            , Element::getBySymbol("C")
            , 3
            , "LEU"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("LEU", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(371)
            , Element::getBySymbol("H")
            , 1
            , "LEU"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("LEU", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(372)
            , Element::getBySymbol("O")
            , 1
            , "LEU"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("LEU", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(373)
            , Element::getBySymbol("H")
            , 1
            , "LEU"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ILE", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(374)
            , Element::getBySymbol("N")
            , 4
            , "ILE"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ILE", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(375)
            , Element::getBySymbol("C")
            , 4
            , "ILE"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ILE", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(376)
            , Element::getBySymbol("C")
            , 3
            , "ILE"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ILE", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(377)
            , Element::getBySymbol("H")
            , 1
            , "ILE"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ILE", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(378)
            , Element::getBySymbol("O")
            , 1
            , "ILE"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ILE", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(379)
            , Element::getBySymbol("H")
            , 1
            , "ILE"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("SER", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(380)
            , Element::getBySymbol("N")
            , 4
            , "SER"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("SER", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(381)
            , Element::getBySymbol("C")
            , 4
            , "SER"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("SER", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(382)
            , Element::getBySymbol("C")
            , 3
            , "SER"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("SER", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(383)
            , Element::getBySymbol("H")
            , 1
            , "SER"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("SER", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(384)
            , Element::getBySymbol("O")
            , 1
            , "SER"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("SER", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(385)
            , Element::getBySymbol("H")
            , 1
            , "SER"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("THR", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(386)
            , Element::getBySymbol("N")
            , 4
            , "THR"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("THR", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(387)
            , Element::getBySymbol("C")
            , 4
            , "THR"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("THR", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(388)
            , Element::getBySymbol("C")
            , 3
            , "THR"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("THR", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(389)
            , Element::getBySymbol("H")
            , 1
            , "THR"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("THR", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(390)
            , Element::getBySymbol("O")
            , 1
            , "THR"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("THR", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(391)
            , Element::getBySymbol("H")
            , 1
            , "THR"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("CYS (-SH)", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(392)
            , Element::getBySymbol("N")
            , 4
            , "CYS (-SH)"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("CYS (-SH)", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(393)
            , Element::getBySymbol("C")
            , 4
            , "CYS (-SH)"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("CYS (-SH)", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(394)
            , Element::getBySymbol("C")
            , 3
            , "CYS (-SH)"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("CYS (-SH)", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(395)
            , Element::getBySymbol("H")
            , 1
            , "CYS (-SH)"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("CYS (-SH)", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(396)
            , Element::getBySymbol("O")
            , 1
            , "CYS (-SH)"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("CYS (-SH)", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(397)
            , Element::getBySymbol("H")
            , 1
            , "CYS (-SH)"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("CYS (-SS)", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(398)
            , Element::getBySymbol("N")
            , 4
            , "CYS (-SS)"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("CYS (-SS)", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(399)
            , Element::getBySymbol("C")
            , 4
            , "CYS (-SS)"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("CYS (-SS)", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(400)
            , Element::getBySymbol("C")
            , 3
            , "CYS (-SS)"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("CYS (-SS)", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(401)
            , Element::getBySymbol("H")
            , 1
            , "CYS (-SS)"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("CYS (-SS)", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(402)
            , Element::getBySymbol("O")
            , 1
            , "CYS (-SS)"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("CYS (-SS)", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(403)
            , Element::getBySymbol("H")
            , 1
            , "CYS (-SS)"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("PRO", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(404)
            , Element::getBySymbol("N")
            , 4
            , "PRO"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("PRO", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(405)
            , Element::getBySymbol("C")
            , 4
            , "PRO"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("PRO", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(406)
            , Element::getBySymbol("C")
            , 3
            , "PRO"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("PRO", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(407)
            , Element::getBySymbol("H")
            , 1
            , "PRO"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("PRO", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(408)
            , Element::getBySymbol("O")
            , 1
            , "PRO"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("PRO", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(409)
            , Element::getBySymbol("H")
            , 1
            , "PRO"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("PRO", "CD", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(410)
            , Element::getBySymbol("C")
            , 4
            , "PRO"
            , "CD"
            , Ordinality::Initial
            );

    if (! Biotype::exists("PRO", "HD", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(411)
            , Element::getBySymbol("H")
            , 1
            , "PRO"
            , "HD"
            , Ordinality::Initial
            );

    if (! Biotype::exists("PHE", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(412)
            , Element::getBySymbol("N")
            , 4
            , "PHE"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("PHE", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(413)
            , Element::getBySymbol("C")
            , 4
            , "PHE"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("PHE", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(414)
            , Element::getBySymbol("C")
            , 3
            , "PHE"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("PHE", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(415)
            , Element::getBySymbol("H")
            , 1
            , "PHE"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("PHE", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(416)
            , Element::getBySymbol("O")
            , 1
            , "PHE"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("PHE", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(417)
            , Element::getBySymbol("H")
            , 1
            , "PHE"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("TYR", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(418)
            , Element::getBySymbol("N")
            , 4
            , "TYR"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("TYR", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(419)
            , Element::getBySymbol("C")
            , 4
            , "TYR"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("TYR", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(420)
            , Element::getBySymbol("C")
            , 3
            , "TYR"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("TYR", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(421)
            , Element::getBySymbol("H")
            , 1
            , "TYR"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("TYR", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(422)
            , Element::getBySymbol("O")
            , 1
            , "TYR"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("TYR", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(423)
            , Element::getBySymbol("H")
            , 1
            , "TYR"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("TRP", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(424)
            , Element::getBySymbol("N")
            , 4
            , "TRP"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("TRP", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(425)
            , Element::getBySymbol("C")
            , 4
            , "TRP"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("TRP", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(426)
            , Element::getBySymbol("C")
            , 3
            , "TRP"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("TRP", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(427)
            , Element::getBySymbol("H")
            , 1
            , "TRP"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("TRP", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(428)
            , Element::getBySymbol("O")
            , 1
            , "TRP"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("TRP", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(429)
            , Element::getBySymbol("H")
            , 1
            , "TRP"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (+)", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(430)
            , Element::getBySymbol("N")
            , 4
            , "HIS (+)"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (+)", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(431)
            , Element::getBySymbol("C")
            , 4
            , "HIS (+)"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (+)", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(432)
            , Element::getBySymbol("C")
            , 3
            , "HIS (+)"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (+)", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(433)
            , Element::getBySymbol("H")
            , 1
            , "HIS (+)"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (+)", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(434)
            , Element::getBySymbol("O")
            , 1
            , "HIS (+)"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (+)", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(435)
            , Element::getBySymbol("H")
            , 1
            , "HIS (+)"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (HD)", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(436)
            , Element::getBySymbol("N")
            , 4
            , "HIS (HD)"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (HD)", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(437)
            , Element::getBySymbol("C")
            , 4
            , "HIS (HD)"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (HD)", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(438)
            , Element::getBySymbol("C")
            , 3
            , "HIS (HD)"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (HD)", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(439)
            , Element::getBySymbol("H")
            , 1
            , "HIS (HD)"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (HD)", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(440)
            , Element::getBySymbol("O")
            , 1
            , "HIS (HD)"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (HD)", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(441)
            , Element::getBySymbol("H")
            , 1
            , "HIS (HD)"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (HE)", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(442)
            , Element::getBySymbol("N")
            , 4
            , "HIS (HE)"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (HE)", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(443)
            , Element::getBySymbol("C")
            , 4
            , "HIS (HE)"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (HE)", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(444)
            , Element::getBySymbol("C")
            , 3
            , "HIS (HE)"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (HE)", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(445)
            , Element::getBySymbol("H")
            , 1
            , "HIS (HE)"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (HE)", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(446)
            , Element::getBySymbol("O")
            , 1
            , "HIS (HE)"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("HIS (HE)", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(447)
            , Element::getBySymbol("H")
            , 1
            , "HIS (HE)"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ASP", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(448)
            , Element::getBySymbol("N")
            , 4
            , "ASP"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ASP", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(449)
            , Element::getBySymbol("C")
            , 4
            , "ASP"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ASP", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(450)
            , Element::getBySymbol("C")
            , 3
            , "ASP"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ASP", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(451)
            , Element::getBySymbol("H")
            , 1
            , "ASP"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ASP", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(452)
            , Element::getBySymbol("O")
            , 1
            , "ASP"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ASP", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(453)
            , Element::getBySymbol("H")
            , 1
            , "ASP"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ASN", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(454)
            , Element::getBySymbol("N")
            , 4
            , "ASN"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ASN", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(455)
            , Element::getBySymbol("C")
            , 4
            , "ASN"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ASN", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(456)
            , Element::getBySymbol("C")
            , 3
            , "ASN"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ASN", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(457)
            , Element::getBySymbol("H")
            , 1
            , "ASN"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ASN", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(458)
            , Element::getBySymbol("O")
            , 1
            , "ASN"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ASN", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(459)
            , Element::getBySymbol("H")
            , 1
            , "ASN"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLU", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(460)
            , Element::getBySymbol("N")
            , 4
            , "GLU"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLU", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(461)
            , Element::getBySymbol("C")
            , 4
            , "GLU"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLU", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(462)
            , Element::getBySymbol("C")
            , 3
            , "GLU"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLU", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(463)
            , Element::getBySymbol("H")
            , 1
            , "GLU"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLU", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(464)
            , Element::getBySymbol("O")
            , 1
            , "GLU"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLU", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(465)
            , Element::getBySymbol("H")
            , 1
            , "GLU"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLN", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(466)
            , Element::getBySymbol("N")
            , 4
            , "GLN"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLN", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(467)
            , Element::getBySymbol("C")
            , 4
            , "GLN"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLN", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(468)
            , Element::getBySymbol("C")
            , 3
            , "GLN"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLN", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(469)
            , Element::getBySymbol("H")
            , 1
            , "GLN"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLN", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(470)
            , Element::getBySymbol("O")
            , 1
            , "GLN"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLN", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(471)
            , Element::getBySymbol("H")
            , 1
            , "GLN"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("MET", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(472)
            , Element::getBySymbol("N")
            , 4
            , "MET"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("MET", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(473)
            , Element::getBySymbol("C")
            , 4
            , "MET"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("MET", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(474)
            , Element::getBySymbol("C")
            , 3
            , "MET"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("MET", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(475)
            , Element::getBySymbol("H")
            , 1
            , "MET"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("MET", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(476)
            , Element::getBySymbol("O")
            , 1
            , "MET"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("MET", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(477)
            , Element::getBySymbol("H")
            , 1
            , "MET"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("LYS", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(478)
            , Element::getBySymbol("N")
            , 4
            , "LYS"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("LYS", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(479)
            , Element::getBySymbol("C")
            , 4
            , "LYS"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("LYS", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(480)
            , Element::getBySymbol("C")
            , 3
            , "LYS"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("LYS", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(481)
            , Element::getBySymbol("H")
            , 1
            , "LYS"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("LYS", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(482)
            , Element::getBySymbol("O")
            , 1
            , "LYS"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("LYS", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(483)
            , Element::getBySymbol("H")
            , 1
            , "LYS"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ARG", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(484)
            , Element::getBySymbol("N")
            , 4
            , "ARG"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ARG", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(485)
            , Element::getBySymbol("C")
            , 4
            , "ARG"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ARG", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(486)
            , Element::getBySymbol("C")
            , 3
            , "ARG"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ARG", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(487)
            , Element::getBySymbol("H")
            , 1
            , "ARG"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ARG", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(488)
            , Element::getBySymbol("O")
            , 1
            , "ARG"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ARG", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(489)
            , Element::getBySymbol("H")
            , 1
            , "ARG"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ORN", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(490)
            , Element::getBySymbol("N")
            , 4
            , "ORN"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ORN", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(491)
            , Element::getBySymbol("C")
            , 4
            , "ORN"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ORN", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(492)
            , Element::getBySymbol("C")
            , 3
            , "ORN"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ORN", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(493)
            , Element::getBySymbol("H")
            , 1
            , "ORN"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ORN", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(494)
            , Element::getBySymbol("O")
            , 1
            , "ORN"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("ORN", "HA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(495)
            , Element::getBySymbol("H")
            , 1
            , "ORN"
            , "HA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("AIB", "N", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(496)
            , Element::getBySymbol("N")
            , 4
            , "AIB"
            , "N"
            , Ordinality::Initial
            );

    if (! Biotype::exists("AIB", "CA", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(497)
            , Element::getBySymbol("C")
            , 4
            , "AIB"
            , "CA"
            , Ordinality::Initial
            );

    if (! Biotype::exists("AIB", "C", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(498)
            , Element::getBySymbol("C")
            , 3
            , "AIB"
            , "C"
            , Ordinality::Initial
            );

    if (! Biotype::exists("AIB", "HN", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(499)
            , Element::getBySymbol("H")
            , 1
            , "AIB"
            , "HN"
            , Ordinality::Initial
            );

    if (! Biotype::exists("AIB", "O", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(500)
            , Element::getBySymbol("O")
            , 1
            , "AIB"
            , "O"
            , Ordinality::Initial
            );

    if (! Biotype::exists("GLY", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(501)
            , Element::getBySymbol("N")
            , 3
            , "GLY"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLY", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(502)
            , Element::getBySymbol("C")
            , 4
            , "GLY"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLY", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(503)
            , Element::getBySymbol("C")
            , 3
            , "GLY"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLY", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(504)
            , Element::getBySymbol("H")
            , 1
            , "GLY"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLY", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(505)
            , Element::getBySymbol("O")
            , 1
            , "GLY"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLY", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(506)
            , Element::getBySymbol("H")
            , 1
            , "GLY"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("ALA", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(507)
            , Element::getBySymbol("N")
            , 3
            , "ALA"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("ALA", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(508)
            , Element::getBySymbol("C")
            , 4
            , "ALA"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("ALA", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(509)
            , Element::getBySymbol("C")
            , 3
            , "ALA"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("ALA", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(510)
            , Element::getBySymbol("H")
            , 1
            , "ALA"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("ALA", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(511)
            , Element::getBySymbol("O")
            , 1
            , "ALA"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("ALA", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(512)
            , Element::getBySymbol("H")
            , 1
            , "ALA"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("VAL", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(513)
            , Element::getBySymbol("N")
            , 3
            , "VAL"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("VAL", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(514)
            , Element::getBySymbol("C")
            , 4
            , "VAL"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("VAL", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(515)
            , Element::getBySymbol("C")
            , 3
            , "VAL"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("VAL", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(516)
            , Element::getBySymbol("H")
            , 1
            , "VAL"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("VAL", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(517)
            , Element::getBySymbol("O")
            , 1
            , "VAL"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("VAL", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(518)
            , Element::getBySymbol("H")
            , 1
            , "VAL"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("LEU", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(519)
            , Element::getBySymbol("N")
            , 3
            , "LEU"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("LEU", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(520)
            , Element::getBySymbol("C")
            , 4
            , "LEU"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("LEU", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(521)
            , Element::getBySymbol("C")
            , 3
            , "LEU"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("LEU", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(522)
            , Element::getBySymbol("H")
            , 1
            , "LEU"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("LEU", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(523)
            , Element::getBySymbol("O")
            , 1
            , "LEU"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("LEU", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(524)
            , Element::getBySymbol("H")
            , 1
            , "LEU"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("ILE", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(525)
            , Element::getBySymbol("N")
            , 3
            , "ILE"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("ILE", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(526)
            , Element::getBySymbol("C")
            , 4
            , "ILE"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("ILE", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(527)
            , Element::getBySymbol("C")
            , 3
            , "ILE"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("ILE", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(528)
            , Element::getBySymbol("H")
            , 1
            , "ILE"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("ILE", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(529)
            , Element::getBySymbol("O")
            , 1
            , "ILE"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("ILE", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(530)
            , Element::getBySymbol("H")
            , 1
            , "ILE"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("SER", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(531)
            , Element::getBySymbol("N")
            , 3
            , "SER"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("SER", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(532)
            , Element::getBySymbol("C")
            , 4
            , "SER"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("SER", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(533)
            , Element::getBySymbol("C")
            , 3
            , "SER"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("SER", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(534)
            , Element::getBySymbol("H")
            , 1
            , "SER"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("SER", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(535)
            , Element::getBySymbol("O")
            , 1
            , "SER"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("SER", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(536)
            , Element::getBySymbol("H")
            , 1
            , "SER"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("THR", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(537)
            , Element::getBySymbol("N")
            , 3
            , "THR"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("THR", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(538)
            , Element::getBySymbol("C")
            , 4
            , "THR"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("THR", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(539)
            , Element::getBySymbol("C")
            , 3
            , "THR"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("THR", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(540)
            , Element::getBySymbol("H")
            , 1
            , "THR"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("THR", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(541)
            , Element::getBySymbol("O")
            , 1
            , "THR"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("THR", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(542)
            , Element::getBySymbol("H")
            , 1
            , "THR"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("CYS (-SH)", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(543)
            , Element::getBySymbol("N")
            , 3
            , "CYS (-SH)"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("CYS (-SH)", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(544)
            , Element::getBySymbol("C")
            , 4
            , "CYS (-SH)"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("CYS (-SH)", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(545)
            , Element::getBySymbol("C")
            , 3
            , "CYS (-SH)"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("CYS (-SH)", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(546)
            , Element::getBySymbol("H")
            , 1
            , "CYS (-SH)"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("CYS (-SH)", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(547)
            , Element::getBySymbol("O")
            , 1
            , "CYS (-SH)"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("CYS (-SH)", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(548)
            , Element::getBySymbol("H")
            , 1
            , "CYS (-SH)"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("CYS (-SS)", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(549)
            , Element::getBySymbol("N")
            , 3
            , "CYS (-SS)"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("CYS (-SS)", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(550)
            , Element::getBySymbol("C")
            , 4
            , "CYS (-SS)"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("CYS (-SS)", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(551)
            , Element::getBySymbol("C")
            , 3
            , "CYS (-SS)"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("CYS (-SS)", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(552)
            , Element::getBySymbol("H")
            , 1
            , "CYS (-SS)"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("CYS (-SS)", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(553)
            , Element::getBySymbol("O")
            , 1
            , "CYS (-SS)"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("CYS (-SS)", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(554)
            , Element::getBySymbol("H")
            , 1
            , "CYS (-SS)"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("PRO", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(555)
            , Element::getBySymbol("N")
            , 3
            , "PRO"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("PRO", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(556)
            , Element::getBySymbol("C")
            , 4
            , "PRO"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("PRO", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(557)
            , Element::getBySymbol("C")
            , 3
            , "PRO"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("PRO", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(558)
            , Element::getBySymbol("O")
            , 1
            , "PRO"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("PRO", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(559)
            , Element::getBySymbol("H")
            , 1
            , "PRO"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("PHE", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(560)
            , Element::getBySymbol("N")
            , 3
            , "PHE"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("PHE", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(561)
            , Element::getBySymbol("C")
            , 4
            , "PHE"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("PHE", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(562)
            , Element::getBySymbol("C")
            , 3
            , "PHE"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("PHE", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(563)
            , Element::getBySymbol("H")
            , 1
            , "PHE"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("PHE", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(564)
            , Element::getBySymbol("O")
            , 1
            , "PHE"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("PHE", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(565)
            , Element::getBySymbol("H")
            , 1
            , "PHE"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("TYR", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(566)
            , Element::getBySymbol("N")
            , 3
            , "TYR"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("TYR", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(567)
            , Element::getBySymbol("C")
            , 4
            , "TYR"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("TYR", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(568)
            , Element::getBySymbol("C")
            , 3
            , "TYR"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("TYR", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(569)
            , Element::getBySymbol("H")
            , 1
            , "TYR"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("TYR", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(570)
            , Element::getBySymbol("O")
            , 1
            , "TYR"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("TYR", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(571)
            , Element::getBySymbol("H")
            , 1
            , "TYR"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("TRP", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(572)
            , Element::getBySymbol("N")
            , 3
            , "TRP"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("TRP", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(573)
            , Element::getBySymbol("C")
            , 4
            , "TRP"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("TRP", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(574)
            , Element::getBySymbol("C")
            , 3
            , "TRP"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("TRP", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(575)
            , Element::getBySymbol("H")
            , 1
            , "TRP"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("TRP", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(576)
            , Element::getBySymbol("O")
            , 1
            , "TRP"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("TRP", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(577)
            , Element::getBySymbol("H")
            , 1
            , "TRP"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (+)", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(578)
            , Element::getBySymbol("N")
            , 3
            , "HIS (+)"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (+)", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(579)
            , Element::getBySymbol("C")
            , 4
            , "HIS (+)"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (+)", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(580)
            , Element::getBySymbol("C")
            , 3
            , "HIS (+)"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (+)", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(581)
            , Element::getBySymbol("H")
            , 1
            , "HIS (+)"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (+)", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(582)
            , Element::getBySymbol("O")
            , 1
            , "HIS (+)"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (+)", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(583)
            , Element::getBySymbol("H")
            , 1
            , "HIS (+)"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (HD)", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(584)
            , Element::getBySymbol("N")
            , 3
            , "HIS (HD)"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (HD)", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(585)
            , Element::getBySymbol("C")
            , 4
            , "HIS (HD)"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (HD)", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(586)
            , Element::getBySymbol("C")
            , 3
            , "HIS (HD)"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (HD)", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(587)
            , Element::getBySymbol("H")
            , 1
            , "HIS (HD)"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (HD)", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(588)
            , Element::getBySymbol("O")
            , 1
            , "HIS (HD)"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (HD)", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(589)
            , Element::getBySymbol("H")
            , 1
            , "HIS (HD)"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (HE)", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(590)
            , Element::getBySymbol("N")
            , 3
            , "HIS (HE)"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (HE)", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(591)
            , Element::getBySymbol("C")
            , 4
            , "HIS (HE)"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (HE)", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(592)
            , Element::getBySymbol("C")
            , 3
            , "HIS (HE)"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (HE)", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(593)
            , Element::getBySymbol("H")
            , 1
            , "HIS (HE)"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (HE)", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(594)
            , Element::getBySymbol("O")
            , 1
            , "HIS (HE)"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("HIS (HE)", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(595)
            , Element::getBySymbol("H")
            , 1
            , "HIS (HE)"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("ASP", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(596)
            , Element::getBySymbol("N")
            , 3
            , "ASP"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("ASP", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(597)
            , Element::getBySymbol("C")
            , 4
            , "ASP"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("ASP", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(598)
            , Element::getBySymbol("C")
            , 3
            , "ASP"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("ASP", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(599)
            , Element::getBySymbol("H")
            , 1
            , "ASP"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("ASP", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(600)
            , Element::getBySymbol("O")
            , 1
            , "ASP"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("ASP", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(601)
            , Element::getBySymbol("H")
            , 1
            , "ASP"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("ASN", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(602)
            , Element::getBySymbol("N")
            , 3
            , "ASN"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("ASN", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(603)
            , Element::getBySymbol("C")
            , 4
            , "ASN"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("ASN", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(604)
            , Element::getBySymbol("C")
            , 3
            , "ASN"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("ASN", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(605)
            , Element::getBySymbol("H")
            , 1
            , "ASN"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("ASN", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(606)
            , Element::getBySymbol("O")
            , 1
            , "ASN"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("ASN", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(607)
            , Element::getBySymbol("H")
            , 1
            , "ASN"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLU", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(608)
            , Element::getBySymbol("N")
            , 3
            , "GLU"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLU", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(609)
            , Element::getBySymbol("C")
            , 4
            , "GLU"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLU", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(610)
            , Element::getBySymbol("C")
            , 3
            , "GLU"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLU", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(611)
            , Element::getBySymbol("H")
            , 1
            , "GLU"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLU", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(612)
            , Element::getBySymbol("O")
            , 1
            , "GLU"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLU", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(613)
            , Element::getBySymbol("H")
            , 1
            , "GLU"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLN", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(614)
            , Element::getBySymbol("N")
            , 3
            , "GLN"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLN", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(615)
            , Element::getBySymbol("C")
            , 4
            , "GLN"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLN", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(616)
            , Element::getBySymbol("C")
            , 3
            , "GLN"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLN", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(617)
            , Element::getBySymbol("H")
            , 1
            , "GLN"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLN", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(618)
            , Element::getBySymbol("O")
            , 1
            , "GLN"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("GLN", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(619)
            , Element::getBySymbol("H")
            , 1
            , "GLN"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("MET", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(620)
            , Element::getBySymbol("N")
            , 3
            , "MET"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("MET", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(621)
            , Element::getBySymbol("C")
            , 4
            , "MET"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("MET", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(622)
            , Element::getBySymbol("C")
            , 3
            , "MET"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("MET", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(623)
            , Element::getBySymbol("H")
            , 1
            , "MET"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("MET", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(624)
            , Element::getBySymbol("O")
            , 1
            , "MET"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("MET", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(625)
            , Element::getBySymbol("H")
            , 1
            , "MET"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("LYS", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(626)
            , Element::getBySymbol("N")
            , 3
            , "LYS"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("LYS", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(627)
            , Element::getBySymbol("C")
            , 4
            , "LYS"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("LYS", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(628)
            , Element::getBySymbol("C")
            , 3
            , "LYS"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("LYS", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(629)
            , Element::getBySymbol("H")
            , 1
            , "LYS"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("LYS", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(630)
            , Element::getBySymbol("O")
            , 1
            , "LYS"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("LYS", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(631)
            , Element::getBySymbol("H")
            , 1
            , "LYS"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("ARG", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(632)
            , Element::getBySymbol("N")
            , 3
            , "ARG"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("ARG", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(633)
            , Element::getBySymbol("C")
            , 4
            , "ARG"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("ARG", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(634)
            , Element::getBySymbol("C")
            , 3
            , "ARG"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("ARG", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(635)
            , Element::getBySymbol("H")
            , 1
            , "ARG"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("ARG", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(636)
            , Element::getBySymbol("O")
            , 1
            , "ARG"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("ARG", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(637)
            , Element::getBySymbol("H")
            , 1
            , "ARG"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("ORN", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(638)
            , Element::getBySymbol("N")
            , 3
            , "ORN"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("ORN", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(639)
            , Element::getBySymbol("C")
            , 4
            , "ORN"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("ORN", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(640)
            , Element::getBySymbol("C")
            , 3
            , "ORN"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("ORN", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(641)
            , Element::getBySymbol("H")
            , 1
            , "ORN"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("ORN", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(642)
            , Element::getBySymbol("O")
            , 1
            , "ORN"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("ORN", "HA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(643)
            , Element::getBySymbol("H")
            , 1
            , "ORN"
            , "HA"
            , Ordinality::Final
            );

    if (! Biotype::exists("AIB", "N", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(644)
            , Element::getBySymbol("N")
            , 3
            , "AIB"
            , "N"
            , Ordinality::Final
            );

    if (! Biotype::exists("AIB", "CA", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(645)
            , Element::getBySymbol("C")
            , 4
            , "AIB"
            , "CA"
            , Ordinality::Final
            );

    if (! Biotype::exists("AIB", "C", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(646)
            , Element::getBySymbol("C")
            , 3
            , "AIB"
            , "C"
            , Ordinality::Final
            );

    if (! Biotype::exists("AIB", "HN", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(647)
            , Element::getBySymbol("H")
            , 1
            , "AIB"
            , "HN"
            , Ordinality::Final
            );

    if (! Biotype::exists("AIB", "OXT", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(648)
            , Element::getBySymbol("O")
            , 1
            , "AIB"
            , "OXT"
            , Ordinality::Final
            );

    if (! Biotype::exists("Adenosine", "O5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1001)
            , Element::getBySymbol("O")
            , 2
            , "Adenosine"
            , "O5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "C5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1002)
            , Element::getBySymbol("C")
            , 4
            , "Adenosine"
            , "C5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "H5'1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1003)
            , Element::getBySymbol("H")
            , 1
            , "Adenosine"
            , "H5'1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "H5'2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1004)
            , Element::getBySymbol("H")
            , 1
            , "Adenosine"
            , "H5'2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "C4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1005)
            , Element::getBySymbol("C")
            , 4
            , "Adenosine"
            , "C4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "H4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1006)
            , Element::getBySymbol("H")
            , 1
            , "Adenosine"
            , "H4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "O4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1007)
            , Element::getBySymbol("O")
            , 2
            , "Adenosine"
            , "O4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "C1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1008)
            , Element::getBySymbol("C")
            , 4
            , "Adenosine"
            , "C1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "H1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1009)
            , Element::getBySymbol("H")
            , 1
            , "Adenosine"
            , "H1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "C3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1010)
            , Element::getBySymbol("C")
            , 4
            , "Adenosine"
            , "C3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "H3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1011)
            , Element::getBySymbol("H")
            , 1
            , "Adenosine"
            , "H3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "C2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1012)
            , Element::getBySymbol("C")
            , 4
            , "Adenosine"
            , "C2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "H2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1013)
            , Element::getBySymbol("H")
            , 1
            , "Adenosine"
            , "H2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "O2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1014)
            , Element::getBySymbol("O")
            , 2
            , "Adenosine"
            , "O2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "HO'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1015)
            , Element::getBySymbol("H")
            , 1
            , "Adenosine"
            , "HO'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "O3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1016)
            , Element::getBySymbol("O")
            , 2
            , "Adenosine"
            , "O3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "N9", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1017)
            , Element::getBySymbol("N")
            , 3
            , "Adenosine"
            , "N9"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "C4", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1018)
            , Element::getBySymbol("C")
            , 3
            , "Adenosine"
            , "C4"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "C5", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1019)
            , Element::getBySymbol("C")
            , 3
            , "Adenosine"
            , "C5"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "N7", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1020)
            , Element::getBySymbol("N")
            , 2
            , "Adenosine"
            , "N7"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "C8", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1021)
            , Element::getBySymbol("C")
            , 3
            , "Adenosine"
            , "C8"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "N3", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1022)
            , Element::getBySymbol("N")
            , 2
            , "Adenosine"
            , "N3"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "C2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1023)
            , Element::getBySymbol("C")
            , 3
            , "Adenosine"
            , "C2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "N1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1024)
            , Element::getBySymbol("N")
            , 2
            , "Adenosine"
            , "N1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "C6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1025)
            , Element::getBySymbol("C")
            , 3
            , "Adenosine"
            , "C6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "H2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1026)
            , Element::getBySymbol("H")
            , 1
            , "Adenosine"
            , "H2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "N6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1027)
            , Element::getBySymbol("N")
            , 3
            , "Adenosine"
            , "N6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "H61", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1028)
            , Element::getBySymbol("H")
            , 1
            , "Adenosine"
            , "H61"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "H62", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1029)
            , Element::getBySymbol("H")
            , 1
            , "Adenosine"
            , "H62"
            , Ordinality::Any
            );

    if (! Biotype::exists("Adenosine", "H8", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1030)
            , Element::getBySymbol("H")
            , 1
            , "Adenosine"
            , "H8"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "O5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1031)
            , Element::getBySymbol("O")
            , 2
            , "Guanosine"
            , "O5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "C5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1032)
            , Element::getBySymbol("C")
            , 4
            , "Guanosine"
            , "C5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "H5'1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1033)
            , Element::getBySymbol("H")
            , 1
            , "Guanosine"
            , "H5'1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "H5'2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1034)
            , Element::getBySymbol("H")
            , 1
            , "Guanosine"
            , "H5'2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "C4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1035)
            , Element::getBySymbol("C")
            , 4
            , "Guanosine"
            , "C4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "H4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1036)
            , Element::getBySymbol("H")
            , 1
            , "Guanosine"
            , "H4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "O4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1037)
            , Element::getBySymbol("O")
            , 2
            , "Guanosine"
            , "O4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "C1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1038)
            , Element::getBySymbol("C")
            , 4
            , "Guanosine"
            , "C1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "H1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1039)
            , Element::getBySymbol("H")
            , 1
            , "Guanosine"
            , "H1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "C3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1040)
            , Element::getBySymbol("C")
            , 4
            , "Guanosine"
            , "C3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "H3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1041)
            , Element::getBySymbol("H")
            , 1
            , "Guanosine"
            , "H3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "C2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1042)
            , Element::getBySymbol("C")
            , 4
            , "Guanosine"
            , "C2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "H2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1043)
            , Element::getBySymbol("H")
            , 1
            , "Guanosine"
            , "H2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "O2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1044)
            , Element::getBySymbol("O")
            , 2
            , "Guanosine"
            , "O2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "HO'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1045)
            , Element::getBySymbol("H")
            , 1
            , "Guanosine"
            , "HO'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "O3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1046)
            , Element::getBySymbol("O")
            , 2
            , "Guanosine"
            , "O3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "N9", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1047)
            , Element::getBySymbol("N")
            , 3
            , "Guanosine"
            , "N9"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "C4", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1048)
            , Element::getBySymbol("C")
            , 3
            , "Guanosine"
            , "C4"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "C5", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1049)
            , Element::getBySymbol("C")
            , 3
            , "Guanosine"
            , "C5"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "N7", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1050)
            , Element::getBySymbol("N")
            , 2
            , "Guanosine"
            , "N7"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "C8", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1051)
            , Element::getBySymbol("C")
            , 3
            , "Guanosine"
            , "C8"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "N3", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1052)
            , Element::getBySymbol("N")
            , 2
            , "Guanosine"
            , "N3"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "C2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1053)
            , Element::getBySymbol("C")
            , 3
            , "Guanosine"
            , "C2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "N1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1054)
            , Element::getBySymbol("N")
            , 3
            , "Guanosine"
            , "N1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "C6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1055)
            , Element::getBySymbol("C")
            , 3
            , "Guanosine"
            , "C6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "H1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1056)
            , Element::getBySymbol("H")
            , 1
            , "Guanosine"
            , "H1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "N2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1057)
            , Element::getBySymbol("N")
            , 3
            , "Guanosine"
            , "N2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "H21", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1058)
            , Element::getBySymbol("H")
            , 1
            , "Guanosine"
            , "H21"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "H22", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1059)
            , Element::getBySymbol("H")
            , 1
            , "Guanosine"
            , "H22"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "O6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1060)
            , Element::getBySymbol("O")
            , 1
            , "Guanosine"
            , "O6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Guanosine", "H8", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1061)
            , Element::getBySymbol("H")
            , 1
            , "Guanosine"
            , "H8"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "O5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1062)
            , Element::getBySymbol("O")
            , 2
            , "Cytidine"
            , "O5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "C5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1063)
            , Element::getBySymbol("C")
            , 4
            , "Cytidine"
            , "C5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "H5'1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1064)
            , Element::getBySymbol("H")
            , 1
            , "Cytidine"
            , "H5'1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "H5'2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1065)
            , Element::getBySymbol("H")
            , 1
            , "Cytidine"
            , "H5'2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "C4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1066)
            , Element::getBySymbol("C")
            , 4
            , "Cytidine"
            , "C4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "H4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1067)
            , Element::getBySymbol("H")
            , 1
            , "Cytidine"
            , "H4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "O4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1068)
            , Element::getBySymbol("O")
            , 2
            , "Cytidine"
            , "O4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "C1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1069)
            , Element::getBySymbol("C")
            , 4
            , "Cytidine"
            , "C1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "H1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1070)
            , Element::getBySymbol("H")
            , 1
            , "Cytidine"
            , "H1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "C3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1071)
            , Element::getBySymbol("C")
            , 4
            , "Cytidine"
            , "C3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "H3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1072)
            , Element::getBySymbol("H")
            , 1
            , "Cytidine"
            , "H3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "C2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1073)
            , Element::getBySymbol("C")
            , 4
            , "Cytidine"
            , "C2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "H2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1074)
            , Element::getBySymbol("H")
            , 1
            , "Cytidine"
            , "H2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "O2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1075)
            , Element::getBySymbol("O")
            , 2
            , "Cytidine"
            , "O2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "HO'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1076)
            , Element::getBySymbol("H")
            , 1
            , "Cytidine"
            , "HO'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "O3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1077)
            , Element::getBySymbol("O")
            , 2
            , "Cytidine"
            , "O3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "N1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1078)
            , Element::getBySymbol("N")
            , 3
            , "Cytidine"
            , "N1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "C2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1079)
            , Element::getBySymbol("C")
            , 3
            , "Cytidine"
            , "C2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "N3", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1080)
            , Element::getBySymbol("N")
            , 2
            , "Cytidine"
            , "N3"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "C4", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1081)
            , Element::getBySymbol("C")
            , 3
            , "Cytidine"
            , "C4"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "C5", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1082)
            , Element::getBySymbol("C")
            , 3
            , "Cytidine"
            , "C5"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "C6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1083)
            , Element::getBySymbol("C")
            , 3
            , "Cytidine"
            , "C6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "O2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1084)
            , Element::getBySymbol("O")
            , 1
            , "Cytidine"
            , "O2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "N4", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1085)
            , Element::getBySymbol("N")
            , 3
            , "Cytidine"
            , "N4"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "H41", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1086)
            , Element::getBySymbol("H")
            , 1
            , "Cytidine"
            , "H41"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "H42", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1087)
            , Element::getBySymbol("H")
            , 1
            , "Cytidine"
            , "H42"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "H5", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1088)
            , Element::getBySymbol("H")
            , 1
            , "Cytidine"
            , "H5"
            , Ordinality::Any
            );

    if (! Biotype::exists("Cytidine", "H6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1089)
            , Element::getBySymbol("H")
            , 1
            , "Cytidine"
            , "H6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "O5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1090)
            , Element::getBySymbol("O")
            , 2
            , "Uridine"
            , "O5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "C5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1091)
            , Element::getBySymbol("C")
            , 4
            , "Uridine"
            , "C5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "H5'1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1092)
            , Element::getBySymbol("H")
            , 1
            , "Uridine"
            , "H5'1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "H5'2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1093)
            , Element::getBySymbol("H")
            , 1
            , "Uridine"
            , "H5'2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "C4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1094)
            , Element::getBySymbol("C")
            , 4
            , "Uridine"
            , "C4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "H4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1095)
            , Element::getBySymbol("H")
            , 1
            , "Uridine"
            , "H4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "O4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1096)
            , Element::getBySymbol("O")
            , 2
            , "Uridine"
            , "O4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "C1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1097)
            , Element::getBySymbol("C")
            , 4
            , "Uridine"
            , "C1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "H1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1098)
            , Element::getBySymbol("H")
            , 1
            , "Uridine"
            , "H1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "C3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1099)
            , Element::getBySymbol("C")
            , 4
            , "Uridine"
            , "C3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "H3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1100)
            , Element::getBySymbol("H")
            , 1
            , "Uridine"
            , "H3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "C2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1101)
            , Element::getBySymbol("C")
            , 4
            , "Uridine"
            , "C2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "H2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1102)
            , Element::getBySymbol("H")
            , 1
            , "Uridine"
            , "H2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "O2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1103)
            , Element::getBySymbol("O")
            , 2
            , "Uridine"
            , "O2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "HO'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1104)
            , Element::getBySymbol("H")
            , 1
            , "Uridine"
            , "HO'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "O3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1105)
            , Element::getBySymbol("O")
            , 2
            , "Uridine"
            , "O3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "N1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1106)
            , Element::getBySymbol("N")
            , 3
            , "Uridine"
            , "N1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "C2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1107)
            , Element::getBySymbol("C")
            , 3
            , "Uridine"
            , "C2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "N3", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1108)
            , Element::getBySymbol("N")
            , 3
            , "Uridine"
            , "N3"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "C4", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1109)
            , Element::getBySymbol("C")
            , 3
            , "Uridine"
            , "C4"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "C5", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1110)
            , Element::getBySymbol("C")
            , 3
            , "Uridine"
            , "C5"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "C6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1111)
            , Element::getBySymbol("C")
            , 3
            , "Uridine"
            , "C6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "O2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1112)
            , Element::getBySymbol("O")
            , 1
            , "Uridine"
            , "O2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "H3", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1113)
            , Element::getBySymbol("H")
            , 1
            , "Uridine"
            , "H3"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "O4", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1114)
            , Element::getBySymbol("O")
            , 1
            , "Uridine"
            , "O4"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "H5", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1115)
            , Element::getBySymbol("H")
            , 1
            , "Uridine"
            , "H5"
            , Ordinality::Any
            );

    if (! Biotype::exists("Uridine", "H6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1116)
            , Element::getBySymbol("H")
            , 1
            , "Uridine"
            , "H6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "O5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1117)
            , Element::getBySymbol("O")
            , 2
            , "Deoxyadenosine"
            , "O5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "C5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1118)
            , Element::getBySymbol("C")
            , 4
            , "Deoxyadenosine"
            , "C5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "H5'1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1119)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyadenosine"
            , "H5'1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "H5'2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1120)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyadenosine"
            , "H5'2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "C4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1121)
            , Element::getBySymbol("C")
            , 4
            , "Deoxyadenosine"
            , "C4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "H4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1122)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyadenosine"
            , "H4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "O4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1123)
            , Element::getBySymbol("O")
            , 2
            , "Deoxyadenosine"
            , "O4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "C1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1124)
            , Element::getBySymbol("C")
            , 4
            , "Deoxyadenosine"
            , "C1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "H1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1125)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyadenosine"
            , "H1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "C3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1126)
            , Element::getBySymbol("C")
            , 4
            , "Deoxyadenosine"
            , "C3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "H3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1127)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyadenosine"
            , "H3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "C2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1128)
            , Element::getBySymbol("C")
            , 4
            , "Deoxyadenosine"
            , "C2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "H2'1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1129)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyadenosine"
            , "H2'1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "H2'2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1130)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyadenosine"
            , "H2'2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "O3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1131)
            , Element::getBySymbol("O")
            , 2
            , "Deoxyadenosine"
            , "O3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "N9", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1132)
            , Element::getBySymbol("N")
            , 3
            , "Deoxyadenosine"
            , "N9"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "C4", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1133)
            , Element::getBySymbol("C")
            , 3
            , "Deoxyadenosine"
            , "C4"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "C5", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1134)
            , Element::getBySymbol("C")
            , 3
            , "Deoxyadenosine"
            , "C5"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "N7", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1135)
            , Element::getBySymbol("N")
            , 2
            , "Deoxyadenosine"
            , "N7"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "C8", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1136)
            , Element::getBySymbol("C")
            , 3
            , "Deoxyadenosine"
            , "C8"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "N3", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1137)
            , Element::getBySymbol("N")
            , 2
            , "Deoxyadenosine"
            , "N3"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "C2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1138)
            , Element::getBySymbol("C")
            , 3
            , "Deoxyadenosine"
            , "C2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "N1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1139)
            , Element::getBySymbol("N")
            , 2
            , "Deoxyadenosine"
            , "N1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "C6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1140)
            , Element::getBySymbol("C")
            , 3
            , "Deoxyadenosine"
            , "C6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "H2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1141)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyadenosine"
            , "H2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "N6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1142)
            , Element::getBySymbol("N")
            , 3
            , "Deoxyadenosine"
            , "N6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "H61", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1143)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyadenosine"
            , "H61"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "H62", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1144)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyadenosine"
            , "H62"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyadenosine", "H8", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1145)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyadenosine"
            , "H8"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "O5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1146)
            , Element::getBySymbol("O")
            , 2
            , "Deoxyguanosine"
            , "O5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "C5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1147)
            , Element::getBySymbol("C")
            , 4
            , "Deoxyguanosine"
            , "C5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "H5'1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1148)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyguanosine"
            , "H5'1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "H5'2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1149)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyguanosine"
            , "H5'2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "C4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1150)
            , Element::getBySymbol("C")
            , 4
            , "Deoxyguanosine"
            , "C4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "H4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1151)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyguanosine"
            , "H4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "O4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1152)
            , Element::getBySymbol("O")
            , 2
            , "Deoxyguanosine"
            , "O4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "C1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1153)
            , Element::getBySymbol("C")
            , 4
            , "Deoxyguanosine"
            , "C1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "H1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1154)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyguanosine"
            , "H1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "C3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1155)
            , Element::getBySymbol("C")
            , 4
            , "Deoxyguanosine"
            , "C3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "H3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1156)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyguanosine"
            , "H3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "C2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1157)
            , Element::getBySymbol("C")
            , 4
            , "Deoxyguanosine"
            , "C2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "H2'1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1158)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyguanosine"
            , "H2'1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "H2'2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1159)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyguanosine"
            , "H2'2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "O3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1160)
            , Element::getBySymbol("O")
            , 2
            , "Deoxyguanosine"
            , "O3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "N9", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1161)
            , Element::getBySymbol("N")
            , 3
            , "Deoxyguanosine"
            , "N9"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "C4", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1162)
            , Element::getBySymbol("C")
            , 3
            , "Deoxyguanosine"
            , "C4"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "C5", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1163)
            , Element::getBySymbol("C")
            , 3
            , "Deoxyguanosine"
            , "C5"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "N7", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1164)
            , Element::getBySymbol("N")
            , 2
            , "Deoxyguanosine"
            , "N7"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "C8", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1165)
            , Element::getBySymbol("C")
            , 3
            , "Deoxyguanosine"
            , "C8"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "N3", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1166)
            , Element::getBySymbol("N")
            , 2
            , "Deoxyguanosine"
            , "N3"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "C2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1167)
            , Element::getBySymbol("C")
            , 3
            , "Deoxyguanosine"
            , "C2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "N1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1168)
            , Element::getBySymbol("N")
            , 3
            , "Deoxyguanosine"
            , "N1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "C6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1169)
            , Element::getBySymbol("C")
            , 3
            , "Deoxyguanosine"
            , "C6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "H1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1170)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyguanosine"
            , "H1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "N2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1171)
            , Element::getBySymbol("N")
            , 3
            , "Deoxyguanosine"
            , "N2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "H21", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1172)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyguanosine"
            , "H21"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "H22", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1173)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyguanosine"
            , "H22"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "O6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1174)
            , Element::getBySymbol("O")
            , 1
            , "Deoxyguanosine"
            , "O6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxyguanosine", "H8", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1175)
            , Element::getBySymbol("H")
            , 1
            , "Deoxyguanosine"
            , "H8"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "O5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1176)
            , Element::getBySymbol("O")
            , 2
            , "Deoxycytidine"
            , "O5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "C5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1177)
            , Element::getBySymbol("C")
            , 4
            , "Deoxycytidine"
            , "C5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "H5'1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1178)
            , Element::getBySymbol("H")
            , 1
            , "Deoxycytidine"
            , "H5'1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "H5'2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1179)
            , Element::getBySymbol("H")
            , 1
            , "Deoxycytidine"
            , "H5'2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "C4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1180)
            , Element::getBySymbol("C")
            , 4
            , "Deoxycytidine"
            , "C4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "H4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1181)
            , Element::getBySymbol("H")
            , 1
            , "Deoxycytidine"
            , "H4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "O4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1182)
            , Element::getBySymbol("O")
            , 2
            , "Deoxycytidine"
            , "O4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "C1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1183)
            , Element::getBySymbol("C")
            , 4
            , "Deoxycytidine"
            , "C1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "H1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1184)
            , Element::getBySymbol("H")
            , 1
            , "Deoxycytidine"
            , "H1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "C3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1185)
            , Element::getBySymbol("C")
            , 4
            , "Deoxycytidine"
            , "C3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "H3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1186)
            , Element::getBySymbol("H")
            , 1
            , "Deoxycytidine"
            , "H3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "C2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1187)
            , Element::getBySymbol("C")
            , 4
            , "Deoxycytidine"
            , "C2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "H2'1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1188)
            , Element::getBySymbol("H")
            , 1
            , "Deoxycytidine"
            , "H2'1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "H2'2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1189)
            , Element::getBySymbol("H")
            , 1
            , "Deoxycytidine"
            , "H2'2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "O3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1190)
            , Element::getBySymbol("O")
            , 2
            , "Deoxycytidine"
            , "O3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "N1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1191)
            , Element::getBySymbol("N")
            , 3
            , "Deoxycytidine"
            , "N1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "C2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1192)
            , Element::getBySymbol("C")
            , 3
            , "Deoxycytidine"
            , "C2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "N3", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1193)
            , Element::getBySymbol("N")
            , 2
            , "Deoxycytidine"
            , "N3"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "C4", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1194)
            , Element::getBySymbol("C")
            , 3
            , "Deoxycytidine"
            , "C4"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "C5", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1195)
            , Element::getBySymbol("C")
            , 3
            , "Deoxycytidine"
            , "C5"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "C6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1196)
            , Element::getBySymbol("C")
            , 3
            , "Deoxycytidine"
            , "C6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "O2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1197)
            , Element::getBySymbol("O")
            , 1
            , "Deoxycytidine"
            , "O2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "N4", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1198)
            , Element::getBySymbol("N")
            , 3
            , "Deoxycytidine"
            , "N4"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "H41", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1199)
            , Element::getBySymbol("H")
            , 1
            , "Deoxycytidine"
            , "H41"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "H42", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1200)
            , Element::getBySymbol("H")
            , 1
            , "Deoxycytidine"
            , "H42"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "H5", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1201)
            , Element::getBySymbol("H")
            , 1
            , "Deoxycytidine"
            , "H5"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxycytidine", "H6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1202)
            , Element::getBySymbol("H")
            , 1
            , "Deoxycytidine"
            , "H6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "O5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1203)
            , Element::getBySymbol("O")
            , 2
            , "Deoxythymidine"
            , "O5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "C5'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1204)
            , Element::getBySymbol("C")
            , 4
            , "Deoxythymidine"
            , "C5'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "H5'1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1205)
            , Element::getBySymbol("H")
            , 1
            , "Deoxythymidine"
            , "H5'1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "H5'2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1206)
            , Element::getBySymbol("H")
            , 1
            , "Deoxythymidine"
            , "H5'2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "C4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1207)
            , Element::getBySymbol("C")
            , 4
            , "Deoxythymidine"
            , "C4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "H4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1208)
            , Element::getBySymbol("H")
            , 1
            , "Deoxythymidine"
            , "H4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "O4'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1209)
            , Element::getBySymbol("O")
            , 2
            , "Deoxythymidine"
            , "O4'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "C1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1210)
            , Element::getBySymbol("C")
            , 4
            , "Deoxythymidine"
            , "C1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "H1'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1211)
            , Element::getBySymbol("H")
            , 1
            , "Deoxythymidine"
            , "H1'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "C3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1212)
            , Element::getBySymbol("C")
            , 4
            , "Deoxythymidine"
            , "C3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "H3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1213)
            , Element::getBySymbol("H")
            , 1
            , "Deoxythymidine"
            , "H3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "C2'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1214)
            , Element::getBySymbol("C")
            , 4
            , "Deoxythymidine"
            , "C2'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "H2'1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1215)
            , Element::getBySymbol("H")
            , 1
            , "Deoxythymidine"
            , "H2'1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "H2'2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1216)
            , Element::getBySymbol("H")
            , 1
            , "Deoxythymidine"
            , "H2'2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "O3'", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1217)
            , Element::getBySymbol("O")
            , 2
            , "Deoxythymidine"
            , "O3'"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "N1", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1218)
            , Element::getBySymbol("N")
            , 3
            , "Deoxythymidine"
            , "N1"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "C2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1219)
            , Element::getBySymbol("C")
            , 3
            , "Deoxythymidine"
            , "C2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "N3", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1220)
            , Element::getBySymbol("N")
            , 3
            , "Deoxythymidine"
            , "N3"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "C4", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1221)
            , Element::getBySymbol("C")
            , 3
            , "Deoxythymidine"
            , "C4"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "C5", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1222)
            , Element::getBySymbol("C")
            , 3
            , "Deoxythymidine"
            , "C5"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "C6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1223)
            , Element::getBySymbol("C")
            , 3
            , "Deoxythymidine"
            , "C6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "O2", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1224)
            , Element::getBySymbol("O")
            , 1
            , "Deoxythymidine"
            , "O2"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "H3", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1225)
            , Element::getBySymbol("H")
            , 1
            , "Deoxythymidine"
            , "H3"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "O4", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1226)
            , Element::getBySymbol("O")
            , 1
            , "Deoxythymidine"
            , "O4"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "C7", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1227)
            , Element::getBySymbol("C")
            , 4
            , "Deoxythymidine"
            , "C7"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "H7", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1228)
            , Element::getBySymbol("H")
            , 1
            , "Deoxythymidine"
            , "H7"
            , Ordinality::Any
            );

    if (! Biotype::exists("Deoxythymidine", "H6", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1229)
            , Element::getBySymbol("H")
            , 1
            , "Deoxythymidine"
            , "H6"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phosphodiester, RNA", "P", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1230)
            , Element::getBySymbol("P")
            , 4
            , "Phosphodiester, RNA"
            , "P"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phosphodiester, RNA", "OP", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1231)
            , Element::getBySymbol("O")
            , 1
            , "Phosphodiester, RNA"
            , "OP"
            , Ordinality::Any
            );

    if (! Biotype::exists("Hydroxyl, RNA", "O5'", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1232)
            , Element::getBySymbol("O")
            , 2
            , "Hydroxyl, RNA"
            , "O5'"
            , Ordinality::Initial
            );

    if (! Biotype::exists("Hydroxyl, RNA", "H5T", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1233)
            , Element::getBySymbol("H")
            , 1
            , "Hydroxyl, RNA"
            , "H5T"
            , Ordinality::Initial
            );

    if (! Biotype::exists("Phosphate, RNA", "O5'", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1234)
            , Element::getBySymbol("O")
            , 2
            , "Phosphate, RNA"
            , "O5'"
            , Ordinality::Initial
            );

    if (! Biotype::exists("Phosphate, RNA", "P", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1235)
            , Element::getBySymbol("P")
            , 4
            , "Phosphate, RNA"
            , "P"
            , Ordinality::Initial
            );

    if (! Biotype::exists("Phosphate, RNA", "OP", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1236)
            , Element::getBySymbol("O")
            , 1
            , "Phosphate, RNA"
            , "OP"
            , Ordinality::Initial
            );

    if (! Biotype::exists("Hydroxyl, RNA", "O3'", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1237)
            , Element::getBySymbol("O")
            , 2
            , "Hydroxyl, RNA"
            , "O3'"
            , Ordinality::Final
            );

    if (! Biotype::exists("Hydroxyl, RNA", "H3T", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1238)
            , Element::getBySymbol("H")
            , 1
            , "Hydroxyl, RNA"
            , "H3T"
            , Ordinality::Final
            );

    if (! Biotype::exists("Phosphate, RNA", "O3'", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1239)
            , Element::getBySymbol("O")
            , 2
            , "Phosphate, RNA"
            , "O3'"
            , Ordinality::Final
            );

    if (! Biotype::exists("Phosphate, RNA", "P", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1240)
            , Element::getBySymbol("P")
            , 4
            , "Phosphate, RNA"
            , "P"
            , Ordinality::Final
            );

    if (! Biotype::exists("Phosphate, RNA", "OP", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1241)
            , Element::getBySymbol("O")
            , 1
            , "Phosphate, RNA"
            , "OP"
            , Ordinality::Final
            );

    if (! Biotype::exists("Phosphodiester, DNA", "P", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1242)
            , Element::getBySymbol("P")
            , 4
            , "Phosphodiester, DNA"
            , "P"
            , Ordinality::Any
            );

    if (! Biotype::exists("Phosphodiester, DNA", "OP", Ordinality::Any) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1243)
            , Element::getBySymbol("O")
            , 1
            , "Phosphodiester, DNA"
            , "OP"
            , Ordinality::Any
            );

    if (! Biotype::exists("Hydroxyl, DNA", "O5'", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1244)
            , Element::getBySymbol("O")
            , 2
            , "Hydroxyl, DNA"
            , "O5'"
            , Ordinality::Initial
            );

    if (! Biotype::exists("Hydroxyl, DNA", "H5T", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1245)
            , Element::getBySymbol("H")
            , 1
            , "Hydroxyl, DNA"
            , "H5T"
            , Ordinality::Initial
            );

    if (! Biotype::exists("Phosphate, DNA", "O5'", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1246)
            , Element::getBySymbol("O")
            , 2
            , "Phosphate, DNA"
            , "O5'"
            , Ordinality::Initial
            );

    if (! Biotype::exists("Phosphate, DNA", "P", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1247)
            , Element::getBySymbol("P")
            , 4
            , "Phosphate, DNA"
            , "P"
            , Ordinality::Initial
            );

    if (! Biotype::exists("Phosphate, DNA", "OP", Ordinality::Initial) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1248)
            , Element::getBySymbol("O")
            , 1
            , "Phosphate, DNA"
            , "OP"
            , Ordinality::Initial
            );

    if (! Biotype::exists("Hydroxyl, DNA", "O3'", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1249)
            , Element::getBySymbol("O")
            , 2
            , "Hydroxyl, DNA"
            , "O3'"
            , Ordinality::Final
            );

    if (! Biotype::exists("Hydroxyl, DNA", "H3T", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1250)
            , Element::getBySymbol("H")
            , 1
            , "Hydroxyl, DNA"
            , "H3T"
            , Ordinality::Final
            );

    if (! Biotype::exists("Phosphate, DNA", "O3'", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1251)
            , Element::getBySymbol("O")
            , 2
            , "Phosphate, DNA"
            , "O3'"
            , Ordinality::Final
            );

    if (! Biotype::exists("Phosphate, DNA", "P", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1252)
            , Element::getBySymbol("P")
            , 4
            , "Phosphate, DNA"
            , "P"
            , Ordinality::Final
            );

    if (! Biotype::exists("Phosphate, DNA", "OP", Ordinality::Final) )
        Biotype::defineTinkerBiotype(
            TinkerBiotypeIndex(1253)
            , Element::getBySymbol("O")
            , 1
            , "Phosphate, DNA"
            , "OP"
            , Ordinality::Final
            );


    // Ions don't have tinker biotype indices

    if (! Biotype::exists("Lithium Ion", "Li+", Ordinality::Any) )
        Biotype::defineBiotype(
            Element::getBySymbol("Li"),
            0, // valence
            "Lithium Ion", "Li+");

    if (! Biotype::exists("Sodium Ion", "Na+", Ordinality::Any) )
        Biotype::defineBiotype(
            Element::getBySymbol("Na"),
            0, // valence
            "Sodium Ion", "Na+");

    if (! Biotype::exists("Potassium Ion", "K+", Ordinality::Any) )
        Biotype::defineBiotype(
            Element::getBySymbol("K"),
            0, // valence
            "Potassium Ion", "K+");

    if (! Biotype::exists("Rubidium Ion", "Rb+", Ordinality::Any) )
        Biotype::defineBiotype(
            Element::getBySymbol("Rb"),
            0, // valence
            "Rubidium Ion", "Rb+");

    if (! Biotype::exists("Cesium Ion", "Cs+", Ordinality::Any) )
        Biotype::defineBiotype(
            Element::getBySymbol("Cs"),
            0, // valence
            "Cesium Ion", "Cs+");

    if (! Biotype::exists("Magnesium Ion", "Mg+2", Ordinality::Any) )
        Biotype::defineBiotype(
            Element::getBySymbol("Mg"),
            0, // valence
            "Magnesium Ion", "Mg+2");

    if (! Biotype::exists("Calcium Ion", "Ca+2", Ordinality::Any) )
        Biotype::defineBiotype(
            Element::getBySymbol("Ca"),
            0, // valence
            "Calcium Ion", "Ca+2");

    if (! Biotype::exists("Zinc Ion", "Zn+2", Ordinality::Any) )
        Biotype::defineBiotype(
            Element::getBySymbol("Zn"),
            0, // valence
            "Zinc Ion", "Zn+2");

    if (! Biotype::exists("Chloride Ion", "Cl-", Ordinality::Any) )
         Biotype::defineBiotype(
           Element::getBySymbol("Cl"),
            0, // valence
           "Chloride Ion", "Cl-");


    popularBiotypesAreInitialized = true;
}


} // namespace SimTK
