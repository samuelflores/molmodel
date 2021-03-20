/* vim: set expandtab ts=4 sts=4 sw=4: */

#ifndef REPRESENTATION_H
#define REPRESENTATION_H

#include <stdexcept>
#include <string>
#include <vector>

namespace Repr {

enum class AtomType : uint8_t {
    UNKNOWN,
    CARBON,
    HYDROGEN,
    NITROGEN,
    OXYGEN,
    SULFUR,
    PHOSPHORUS,
    IRON,
    ALL
};

enum class ResidueProp : uint8_t {
    UNKNOWN,
    ACIDIC,
    BASIC,
    HYDROPHOBIC,
    POLAR
};

enum class ResidueType : uint8_t {
    UNKNOWN,
    /* Aminoacids */
    ALANINE,
    ARGININE,
    ASPARAGINE,
    ASPARTIC_ACID,
    CYSTEINE,
    GLUTAMIC_ACID,
    GLUTAMINE,
    GLYCINE,
    HISTIDINE,
    ISOLEUCINE,
    LEUCINE,
    LYSINE,
    METHIONINE,
    PHENYLALANINE,
    PROLINE,
    SERINE,
    THREONINE,
    TRYPTOPHAN,
    TYROSINE,
    VALINE,
    DISULPHIDEBRIDGEDCYSTEINE,

    /* RNA nucleotides */
    ADENOSINE,
    GUANOSINE,
    CYTOSINE,
    URIDINE,
    THYMINE,
    URIDINE2,
    CYTOSINE2,
    GUANOSINE2,
    URIDINE3,
    CYTOSINE3,
    GUANOSINE3,
    URIDINE4,
    ADENOSINE2,
    GUANOSINE4,
    GUANOSINE5,

    /* DNA deoxynucleotides */
    DEOXYADENOSINE,
    DEOXYGUANOSINE,
    DEOXYCYTOSINE,
    DEOXYTHYMINE,

    /* Misc */
    SOLV,
    ADP,
    CHLORINE,
    HOH
};

class RepresentationException : std::runtime_error {
    using std::runtime_error::runtime_error;
};

class ResidueSpecifier {
public:
#ifdef _MSC_VER
    ResidueSpecifier();
    ResidueSpecifier(const char longName[28], const char abbrevName[4], const char shortName, const ResidueProp prop, const ResidueType type);
#endif // _MSC_VER

    const char longName[28];
    const char abbrevName[4];
    const char shortName;
    const ResidueProp prop;
    const ResidueType type;

    ResidueSpecifier & operator=(const ResidueSpecifier &other);
};

class AtomName {
public:
    AtomName() :
        type{AtomType::UNKNOWN}
    {
        name[0] = '\0';
    }
    explicit AtomName(const std::string &name);

    char name[11];
    AtomType type;
};

#define CHAIN_ID_LEN 4


class Atom {
private:
    char chain_id[CHAIN_ID_LEN];

public:
    int id;
    int orig_id;
    int res_seq;
    int num;
    float temp;
    double posX;
    double posY;
    double posZ;
    AtomName name;
    char insertion_code;

    ResidueType resType;
    ResidueProp resProp;
    bool het;

    const char * getChainId() const {
        return chain_id;
    }

    void setChainId(const std::string &id);
};

class Residue {
public:
    std::vector<Atom> atoms;
    int id;
    ResidueType type;
    ResidueProp prop;
    char insertion_code;
};

class Chain {
private:
    char id[CHAIN_ID_LEN];

public:
    explicit Chain(const std::string &id);

    std::vector<Residue> residues;

    const char * getId() const {
        return id;
    }
};

class Model {
public:
    std::string name;
    std::vector<Chain> chains;
};

using Structure = std::vector<Model>;

AtomType getAtomType(const std::string &name);
ResidueProp getResidueProp(const ResidueType type);
const ResidueSpecifier & getResidueSpecifier(const ResidueType type);
ResidueType getResidueType(const std::string &name);

inline
bool residueIsDNA(const ResidueType type) {
    const auto n = std::underlying_type_t<ResidueType>(type);
    // For reference, look at the list of residue types in Representation.cpp
    return std::underlying_type_t<ResidueType>(ResidueType::DEOXYADENOSINE) <= n &&
           n <= std::underlying_type_t<ResidueType>(ResidueType::DEOXYTHYMINE);
}

inline
bool residueIsProtein(const ResidueType type) {
    const auto n = std::underlying_type_t<ResidueType>(type);
    return std::underlying_type_t<ResidueType>(ResidueType::ALANINE) <= n &&
           n <= std::underlying_type_t<ResidueType>(ResidueType::DISULPHIDEBRIDGEDCYSTEINE);
}

inline
bool residueIsRNA(const ResidueType type) {
    const auto n = std::underlying_type_t<ResidueType>(type);
    // For reference, look at the list of residue types in Representation.cpp
    return std::underlying_type_t<ResidueType>(ResidueType::ADENOSINE) <= n &&
           n <= std::underlying_type_t<ResidueType>(ResidueType::URIDINE);
}

template <typename CharT>
std::basic_string<CharT> rtToString(const ResidueType type);

template <>
inline
std::basic_string<char> rtToString<char>(const ResidueType type) {
    return std::to_string(std::underlying_type_t<ResidueType>(type));
}

template <>
inline
std::basic_string<wchar_t> rtToString<wchar_t>(const ResidueType type) {
    return std::to_wstring(std::underlying_type_t<ResidueType>(type));
}

template <typename CharT>
inline
std::basic_ostream<CharT, std::char_traits<CharT>> & operator<<(std::basic_ostream<CharT, std::char_traits<CharT>> &stm, const ResidueType &type) {
    stm << rtToString<CharT>(type);
    return stm;
}

} // namespace Repr

#endif // REPRESENTATION_H
