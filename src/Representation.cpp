/* vim: set expandtab ts=4 sts=4 sw=4: */

#include "Representation.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstring>
#include <type_traits>

#define CHK_IDENTIFIER(idf) \
    assert(idf.length() > 0); \
    assert(idf.front() != ' '); \
    assert(idf.back() != ' ')

#define MK_RESIDUE_SPECIFIER(longName, abbrevName, shortName, prop, type) \
    arr.at(std::underlying_type_t<ResidueType>(ResidueType::type)) = ResidueSpecifier{longName, abbrevName, shortName, ResidueProp::prop, ResidueType::type}

template <size_t N>
void charArrayFromString(char dst[N], const std::string &src) {
    if (src.length() < 1) {
        dst[0] = '\0';
        return;
    }

    size_t idx = 0;
    for (; idx < src.length(); idx++) {
        if (idx >= N) {
#ifndef MOLMODEL_REPR_TRUNCATE
            std::string msg = "Input string is too long to fit into array of size " + std::to_string(N);
            throw new Repr::RepresentationException{msg};
#else
            break;
#endif // REPR_TRUNCATE
        }
        dst[idx] = src[idx];
    }

    dst[idx >= N ? N - 1 : idx] = '\0';
}

namespace Repr {

ResidueSpecifier & ResidueSpecifier::operator=(const ResidueSpecifier &other) {
    std::strcpy((char *)this->longName, other.longName);
    std::strcpy((char *)this->abbrevName, other.abbrevName);
    const_cast<char&>(this->shortName) = other.shortName;
    const_cast<ResidueProp&>(this->prop) = other.prop;
    const_cast<ResidueType&>(this->type) = other.type;

    return *this;
}

using ResidueSpecifiers = std::array<ResidueSpecifier, 45>;

ResidueSpecifiers initResidueSpecifiers() {
    ResidueSpecifiers arr{};

    MK_RESIDUE_SPECIFIER("unknown",        "unk", ' ', UNKNOWN, UNKNOWN);
    MK_RESIDUE_SPECIFIER("alanine",        "ala", 'a', HYDROPHOBIC, ALANINE),
    MK_RESIDUE_SPECIFIER("arginine",       "arg", 'r', BASIC, ARGININE),
    MK_RESIDUE_SPECIFIER("asparagine",     "asn", 'n', POLAR, ASPARAGINE),
    MK_RESIDUE_SPECIFIER("aspartic acid",  "asp", 'd', ACIDIC, ASPARTIC_ACID),
    MK_RESIDUE_SPECIFIER("cysteine",       "cys", 'c', HYDROPHOBIC, CYSTEINE),
    MK_RESIDUE_SPECIFIER("glutamic acid",  "glu", 'e', ACIDIC, GLUTAMIC_ACID),
    MK_RESIDUE_SPECIFIER("glutamine",      "gln", 'q', POLAR, GLUTAMINE),
    MK_RESIDUE_SPECIFIER("glycine",        "gly", 'g', HYDROPHOBIC, GLYCINE),
    MK_RESIDUE_SPECIFIER("histidine",      "his", 'h', BASIC, HISTIDINE),
    MK_RESIDUE_SPECIFIER("isoleucine",     "ile", 'i', HYDROPHOBIC, ISOLEUCINE),
    MK_RESIDUE_SPECIFIER("leucine",        "leu", 'l', HYDROPHOBIC, LEUCINE),
    MK_RESIDUE_SPECIFIER("lysine",         "lys", 'k', BASIC, LYSINE),
    MK_RESIDUE_SPECIFIER("methionine",     "met", 'm', HYDROPHOBIC, METHIONINE),
    MK_RESIDUE_SPECIFIER("phenylalanine",  "phe", 'f', HYDROPHOBIC, PHENYLALANINE),
    MK_RESIDUE_SPECIFIER("proline",        "pro", 'p', HYDROPHOBIC, PROLINE),
    MK_RESIDUE_SPECIFIER("serine",         "ser", 's', POLAR, SERINE),
    MK_RESIDUE_SPECIFIER("threonine",      "thr", 't', POLAR, THREONINE),
    MK_RESIDUE_SPECIFIER("tryptophan",     "trp", 'w', HYDROPHOBIC, TRYPTOPHAN),
    MK_RESIDUE_SPECIFIER("tyrosine",       "tyr", 'y', HYDROPHOBIC, TYROSINE),
    MK_RESIDUE_SPECIFIER("valine",         "val", 'v', HYDROPHOBIC, VALINE),
    MK_RESIDUE_SPECIFIER("solvent",        "sol", 'h', UNKNOWN, SOLV),
    MK_RESIDUE_SPECIFIER("adp",            "adp", 'h', UNKNOWN, ADP),
    MK_RESIDUE_SPECIFIER("cl",             "cli", 'h', UNKNOWN, CHLORINE),
    MK_RESIDUE_SPECIFIER("hoh",            "hoh", 'h', UNKNOWN, HOH),
    MK_RESIDUE_SPECIFIER("adenosine",      "adn", 'A', UNKNOWN, ADENOSINE),
    MK_RESIDUE_SPECIFIER("guanosine",      "gua", 'G', UNKNOWN, GUANOSINE),
    MK_RESIDUE_SPECIFIER("cytosine",       "cyt", 'C', UNKNOWN, CYTOSINE),
    MK_RESIDUE_SPECIFIER("uridine",        "ura", 'U', UNKNOWN, URIDINE),
    MK_RESIDUE_SPECIFIER("thymine",        "thy", 'T', UNKNOWN, THYMINE),
    MK_RESIDUE_SPECIFIER("uridine",        "h2u", 'U', UNKNOWN, URIDINE2),
    MK_RESIDUE_SPECIFIER("cytosine",       "omc", 'C', UNKNOWN, CYTOSINE2),
    MK_RESIDUE_SPECIFIER("guanosine",      "omg", 'G', UNKNOWN, GUANOSINE2),
    MK_RESIDUE_SPECIFIER("uridine",        "psu", 'U', UNKNOWN, URIDINE3),
    MK_RESIDUE_SPECIFIER("cytosine",       "5mc", 'C', UNKNOWN, CYTOSINE3),
    MK_RESIDUE_SPECIFIER("guanosine",      "7mg", 'G', UNKNOWN, GUANOSINE3),
    MK_RESIDUE_SPECIFIER("uridine",        "5mu", 'U', UNKNOWN, URIDINE4),
    MK_RESIDUE_SPECIFIER("adenosine",      "1ma", 'A', UNKNOWN, ADENOSINE2),
    MK_RESIDUE_SPECIFIER("guanosine",      "2mg", 'G', UNKNOWN, GUANOSINE4),
    MK_RESIDUE_SPECIFIER("guanosine",      "m2g", 'G', UNKNOWN, GUANOSINE5),

    // scf added DNA residues
    MK_RESIDUE_SPECIFIER("deoxyadenosine", "da",  'A', UNKNOWN, DEOXYADENOSINE),
    MK_RESIDUE_SPECIFIER("deoxyguanosine", "dg",  'G', UNKNOWN, DEOXYGUANOSINE),
    MK_RESIDUE_SPECIFIER("deoxycytosine",  "dc",  'C', UNKNOWN, DEOXYCYTOSINE),
    MK_RESIDUE_SPECIFIER("deoxythymidine", "dt",  'T', UNKNOWN, DEOXYTHYMINE),

    //
    MK_RESIDUE_SPECIFIER("disulphidebridgedcysteine", "cyx", 'x', HYDROPHOBIC, DISULPHIDEBRIDGEDCYSTEINE);

    return arr;
}
static const auto RESIDUE_SPECIFIERS = initResidueSpecifiers();

AtomName::AtomName(const std::string &name) :
    type{getAtomType(name)}
{
    charArrayFromString<sizeof(this->name)>(this->name, name);
}

AtomType getAtomType(const std::string &name) {
    CHK_IDENTIFIER(name);

    if (name.length() > 1)
        return AtomType::UNKNOWN;

    switch (name[0]) {
    case 'H':
        return AtomType::HYDROGEN;
    case 'C':
        return AtomType::CARBON;
    case 'O':
        return AtomType::OXYGEN;
    case 'N':
        return AtomType::NITROGEN;
    case 'P':
        return AtomType::PHOSPHORUS;
    case 'S':
        return AtomType::SULFUR;
    default:
        return AtomType::UNKNOWN;
    }
}

void Atom::setChainId(const std::string &id) {
    charArrayFromString<sizeof(chain_id)>(chain_id, id);
}

Chain::Chain(const std::string &id) {
    charArrayFromString<sizeof(this->id)>(this->id, id);
}

ResidueProp getResidueProp(const ResidueType type) {
    return RESIDUE_SPECIFIERS[std::underlying_type_t<decltype(type)>(type)].prop;
}

const ResidueSpecifier & getResidueSpecifier(const ResidueType type) {
    return RESIDUE_SPECIFIERS[std::underlying_type_t<decltype(type)>(type)];
}

ResidueType getResidueType(const std::string &name) {
    CHK_IDENTIFIER(name);

    std::string lcaseName{};
    std::transform(name.cbegin(), name.cend(), std::back_inserter(lcaseName), ::tolower);

    const auto len = name.length();

    if (len == 1) {
        for (const auto &it : RESIDUE_SPECIFIERS) {
            if (lcaseName[0] == it.shortName)
                return it.type;
        }
    } else if (len <= 3) {
        for (const auto &it : RESIDUE_SPECIFIERS) {
            if (lcaseName.compare(it.abbrevName) == 0)
                return it.type;
        }

        // No match, check for GROMACS modified residue names
        std::string groName = "__c";
        for (const auto &it : RESIDUE_SPECIFIERS) {
            groName[0] = it.longName[0];
            groName[1] = it.longName[2];

            if (lcaseName.compare(groName) == 0)
                return it.type;
        }
    } else {
        for (const auto &it : RESIDUE_SPECIFIERS) {
            if (lcaseName.compare(it.longName) == 0)
                return it.type;
        }
    }

    return ResidueType::UNKNOWN;
}

} // namespace Repr
