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
#define DO_INSTANTIATE_ELEMENT_PIMPL_HANDLE
#include "molmodel/internal/Element.h"
#undef DO_INSTANTIATE_ELEMENT_PIMPL_HANDLE

#include "ElementRep.h"
// #include "SimTKcommon/internal/PrivateImplementation_Defs.h"

#include <algorithm>
#include <array>
#include <map>
#include <mutex>
#include <stdexcept>

using namespace std;


namespace SimTK {

// template class PIMPLImplementation<Element,ElementRep>; // explicit instantiation


///////////////
/// Element ///
///////////////

// template class PIMPLHandle<Element,ElementRep>; // explicit instantiation


#define NUM_ELEMENTS 116
static std::array<Element *, NUM_ELEMENTS> ELEMENTS;
static Element * _Deuterium = new Element::Deuterium();

static std::mutex INIT_MTX;

static Element * makeElement(int atomicNumber) {
    switch (atomicNumber) {
        case 1: return new Element::Hydrogen();
        case 2: return new Element::Helium();
        case 3: return new Element::Lithium();
        case 4: return new Element::Beryllium();
        case 5: return new Element::Boron();
        case 6: return new Element::Carbon();
        case 7: return new Element::Nitrogen();
        case 8: return new Element::Oxygen();
        case 9: return new Element::Fluorine();
        case 10: return new Element::Neon();
        case 11: return new Element::Sodium();
        case 12: return new Element::Magnesium();
        case 13: return new Element::Aluminum();
	case 14: return new Element::Silicon();
        case 15: return new Element::Phosphorus();
	case 16: return new Element::Sulfur();
	case 17: return new Element::Chlorine();
	case 18: return new Element::Argon();
        case 19: return new Element::Potassium();
	case 20: return new Element::Calcium();
	case 21: return new Element::Scandium();
	case 22: return new Element::Titanium();
	case 23: return new Element::Vanadium();
	case 24: return new Element::Chromium();
	case 25: return new Element::Manganese();
	case 26: return new Element::Iron();
	case 27: return new Element::Cobalt();
	case 28: return new Element::Nickel();
	case 29: return new Element::Copper();
	case 30: return new Element::Zinc();
	case 31: return new Element::Gallium();
	case 32: return new Element::Germanium();
	case 33: return new Element::Arsenic();
	case 34: return new Element::Selenium();
	case 35: return new Element::Bromine();
	case 36: return new Element::Krypton();
	case 37: return new Element::Rubidium();
	case 38: return new Element::Strontium();
	case 39: return new Element::Yttrium();
	case 40: return new Element::Zirconium();
	case 41: return new Element::Niobium();
	case 42: return new Element::Molybdenum();
	case 43: return new Element::Technetium();
	case 44: return new Element::Ruthenium();
	case 45: return new Element::Rhodium();
	case 46: return new Element::Palladium();
	case 47: return new Element::Silver();
	case 48: return new Element::Cadmium();
	case 49: return new Element::Indium();
	case 50: return new Element::Tin();
	case 51: return new Element::Antimony();
	case 52: return new Element::Tellurium();
	case 53: return new Element::Iodine();
	case 54: return new Element::Xenon();
	case 55: return new Element::Cesium();
	case 56: return new Element::Barium();
	case 57: return new Element::Lanthanum();
	case 58: return new Element::Cerium();
	case 59: return new Element::Praseodymium();
	case 60: return new Element::Neodymium();
	case 61: return new Element::Promethium();
	case 62: return new Element::Samarium();
	case 63: return new Element::Europium();
	case 64: return new Element::Gadolinium();
	case 65: return new Element::Terbium();
	case 66: return new Element::Dysprosium();
	case 67: return new Element::Holmium();
	case 68: return new Element::Erbium();
	case 69: return new Element::Thulium();
	case 70: return new Element::Ytterbium();
	case 71: return new Element::Lutetium();
	case 72: return new Element::Hafnium();
	case 73: return new Element::Tantalum();
	case 74: return new Element::Tungsten();
	case 75: return new Element::Rhenium();
	case 76: return new Element::Osmium();
	case 77: return new Element::Iridium();
	case 78: return new Element::Platinum();
	case 79: return new Element::Gold();
	case 80: return new Element::Mercury();
	case 81: return new Element::Thallium();
	case 82: return new Element::Lead();
	case 83: return new Element::Bismuth();
	case 84: return new Element::Polonium();
	case 85: return new Element::Astatine();
	case 86: return new Element::Radon();
	case 87: return new Element::Francium();
	case 88: return new Element::Radium();
	case 89: return new Element::Actinium();
	case 90: return new Element::Thorium();
	case 91: return new Element::Protactinium();
	case 92: return new Element::Uranium();
	case 93: return new Element::Neptunium();
	case 94: return new Element::Plutonium();
	case 95: return new Element::Americium();
	case 96: return new Element::Curium();
	case 97: return new Element::Berkelium();
	case 98: return new Element::Californium();
	case 99: return new Element::Einsteinium();
	case 100: return new Element::Fermium();
	case 101: return new Element::Mendelevium();
	case 102: return new Element::Nobelium();
	case 103: return new Element::Lawrencium();
	case 104: return new Element::Rutherfordium();
	case 105: return new Element::Dubnium();
	case 106: return new Element::Seaborgium();
	case 107: return new Element::Bohrium();
	case 108: return new Element::Hassium();
	case 109: return new Element::Meitnerium();
	case 110: return new Element::Darmstadtium();
	case 111: return new Element::Roentgenium();
	case 112: return new Element::Ununbium();
	case 113: return new Element::Ununtrium();
	case 114: return new Element::Ununquadium();
	case 115: return new Element::Ununpentium();
	case 116: return new Element::Ununhexium();
    }

    throw std::out_of_range("Atomic number " + std::to_string(atomicNumber) + " does not correspond to any known element");
}

/* static */ const Element * Element::getByAtomicNumber(int atomicNumber) {
    if (atomicNumber < 1 || atomicNumber >= NUM_ELEMENTS)
        throw std::out_of_range("Atomic number " + std::to_string(atomicNumber) + " does not correspond to any known element");

    size_t idx = atomicNumber - 1;

    std::lock_guard<std::mutex> lk{INIT_MTX};
    if (ELEMENTS[idx] == nullptr)
        ELEMENTS[idx] = makeElement(atomicNumber);

    return ELEMENTS[idx];
}

/* static */ const Element * Element::getByName(Name name) {
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);

    if      (name == "hydrogen")  return getByAtomicNumber(1);
    else if (name == "deuterium")  return _Deuterium;
    else if (name == "helium") return getByAtomicNumber(2);
    else if (name == "lithium") return getByAtomicNumber(3);
    else if (name == "beryllium") return getByAtomicNumber(4);
    else if (name == "boron")  return getByAtomicNumber(5);
    else if (name == "carbon")  return getByAtomicNumber(6);
    else if (name == "nitrogen")  return getByAtomicNumber(7);
    else if (name == "oxygen")  return getByAtomicNumber(8);
    else if (name == "fluorine")  return getByAtomicNumber(9);
    else if (name == "neon") return getByAtomicNumber(10);
    else if (name == "sodium") return getByAtomicNumber(11);
    else if (name == "magnesium") return getByAtomicNumber(12);
    else if (name == "aluminum" || name == "aluminium") return getByAtomicNumber(13);
    else if (name == "silicon") return getByAtomicNumber(14);
    else if (name == "phosphorus")  return getByAtomicNumber(15);
    else if (name == "sulfur")  return getByAtomicNumber(16);
    else if (name == "chlorine") return getByAtomicNumber(17);
    else if (name == "argon") return getByAtomicNumber(18);
    else if (name == "potassium")  return getByAtomicNumber(19);
    else if (name == "calcium") return getByAtomicNumber(20);
    else if (name == "scandium") return getByAtomicNumber(21);
    else if (name == "titanium") return getByAtomicNumber(22);
    else if (name == "vanadium")  return getByAtomicNumber(23);
    else if (name == "chromium") return getByAtomicNumber(24);
    else if (name == "manganese") return getByAtomicNumber(25);
    else if (name == "iron") return getByAtomicNumber(26);
    else if (name == "cobalt") return getByAtomicNumber(27);
    else if (name == "nickel") return getByAtomicNumber(28);
    else if (name == "copper") return getByAtomicNumber(29);
    else if (name == "zinc") return getByAtomicNumber(30);
    else if (name == "gallium") return getByAtomicNumber(31);
    else if (name == "germanium") return getByAtomicNumber(32);
    else if (name == "arsenic") return getByAtomicNumber(33);
    else if (name == "selenium") return getByAtomicNumber(34);
    else if (name == "bromine") return getByAtomicNumber(35);
    else if (name == "krypton") return getByAtomicNumber(36);
    else if (name == "rubidium") return getByAtomicNumber(37);
    else if (name == "strontium") return getByAtomicNumber(38);
    else if (name == "yttrium")  return getByAtomicNumber(39);
    else if (name == "zirconium") return getByAtomicNumber(40);
    else if (name == "niobium") return getByAtomicNumber(41);
    else if (name == "molybdenum") return getByAtomicNumber(42);
    else if (name == "technetium") return getByAtomicNumber(43);
    else if (name == "ruthenium") return getByAtomicNumber(44);
    else if (name == "rhodium") return getByAtomicNumber(45);
    else if (name == "palladium") return getByAtomicNumber(46);
    else if (name == "silver") return getByAtomicNumber(47);
    else if (name == "cadmium") return getByAtomicNumber(48);
    else if (name == "indium") return getByAtomicNumber(49);
    else if (name == "tin") return getByAtomicNumber(50);
    else if (name == "antimony") return getByAtomicNumber(51);
    else if (name == "tellurium") return getByAtomicNumber(52);
    else if (name == "iodine")  return getByAtomicNumber(53);
    else if (name == "xenon") return getByAtomicNumber(54);
    else if (name == "cesium") return getByAtomicNumber(55);
    else if (name == "barium") return getByAtomicNumber(56);
    else if (name == "lanthanum") return getByAtomicNumber(57);
    else if (name == "cerium") return getByAtomicNumber(58);
    else if (name == "praseodymium") return getByAtomicNumber(59);
    else if (name == "neodymium") return getByAtomicNumber(60);
    else if (name == "promethium") return getByAtomicNumber(61);
    else if (name == "samarium") return getByAtomicNumber(62);
    else if (name == "europium") return getByAtomicNumber(63);
    else if (name == "gadolinium") return getByAtomicNumber(64);
    else if (name == "terbium") return getByAtomicNumber(65);
    else if (name == "dysprosium") return getByAtomicNumber(66);
    else if (name == "holmium") return getByAtomicNumber(67);
    else if (name == "erbium") return getByAtomicNumber(68);
    else if (name == "thulium") return getByAtomicNumber(69);
    else if (name == "ytterbium") return getByAtomicNumber(70);
    else if (name == "lutetium") return getByAtomicNumber(71);
    else if (name == "hafnium") return getByAtomicNumber(72);
    else if (name == "tantalum") return getByAtomicNumber(73);
    else if (name == "tungsten")  return getByAtomicNumber(74);
    else if (name == "rhenium") return getByAtomicNumber(75);
    else if (name == "osmium") return getByAtomicNumber(76);
    else if (name == "iridium") return getByAtomicNumber(77);
    else if (name == "platinum") return getByAtomicNumber(78);
    else if (name == "gold") return getByAtomicNumber(79);
    else if (name == "mercury") return getByAtomicNumber(80);
    else if (name == "thallium") return getByAtomicNumber(81);
    else if (name == "lead") return getByAtomicNumber(82);
    else if (name == "bismuth") return getByAtomicNumber(83);
    else if (name == "polonium") return getByAtomicNumber(84);
    else if (name == "astatine") return getByAtomicNumber(85);
    else if (name == "radon") return getByAtomicNumber(86);
    else if (name == "francium") return getByAtomicNumber(87);
    else if (name == "radium") return getByAtomicNumber(88);
    else if (name == "actinium") return getByAtomicNumber(89);
    else if (name == "thorium") return getByAtomicNumber(90);
    else if (name == "protactinium") return getByAtomicNumber(91);
    else if (name == "uranium")  return getByAtomicNumber(92);
    else if (name == "neptunium") return getByAtomicNumber(93);
    else if (name == "plutonium") return getByAtomicNumber(94);
    else if (name == "americium") return getByAtomicNumber(95);
    else if (name == "curium") return getByAtomicNumber(96);
    else if (name == "berkelium") return getByAtomicNumber(97);
    else if (name == "californium") return getByAtomicNumber(98);
    else if (name == "einsteinium") return getByAtomicNumber(99);
    else if (name == "fermium") return getByAtomicNumber(100);
    else if (name == "mendelevium") return getByAtomicNumber(101);
    else if (name == "nobelium") return getByAtomicNumber(102);
    else if (name == "lawrencium") return getByAtomicNumber(103);
    else if (name == "rutherfordium") return getByAtomicNumber(104);
    else if (name == "dubnium") return getByAtomicNumber(105);
    else if (name == "seaborgium") return getByAtomicNumber(106);
    else if (name == "bohrium") return getByAtomicNumber(107);
    else if (name == "hassium") return getByAtomicNumber(108);
    else if (name == "meitnerium") return getByAtomicNumber(109);
    else if (name == "darmstadtium") return getByAtomicNumber(110);
    else if (name == "roentgenium") return getByAtomicNumber(111);
    else if (name == "ununbium" || name == "copernicium") return getByAtomicNumber(112);
    else if (name == "ununtrium" || name == "nihonium") return getByAtomicNumber(113);
    else if (name == "ununquadium" || name == "flerovium") return getByAtomicNumber(114);
    else if (name == "ununpentium" || name == "moscovium") return getByAtomicNumber(115);
    else if (name == "Ununhexium" || name == "livermorium") return getByAtomicNumber(116);

    throw std::out_of_range("Invalid element name " + name);
}

/* static */ const Element * Element::getBySymbol(const SimTK::String& symbol) {
    if      (symbol == "H")  return getByAtomicNumber(1);
    else if (symbol == "D")  return _Deuterium;
    else if (symbol == "He") return getByAtomicNumber(2);
    else if (symbol == "Li") return getByAtomicNumber(3);
    else if (symbol == "Be") return getByAtomicNumber(4);
    else if (symbol == "B")  return getByAtomicNumber(5);
    else if (symbol == "C")  return getByAtomicNumber(6);
    else if (symbol == "N")  return getByAtomicNumber(7);
    else if (symbol == "O")  return getByAtomicNumber(8);
    else if (symbol == "F")  return getByAtomicNumber(9);
    else if (symbol == "Ne") return getByAtomicNumber(10);
    else if (symbol == "Na") return getByAtomicNumber(11);
    else if (symbol == "Mg") return getByAtomicNumber(12);
    else if (symbol == "Al") return getByAtomicNumber(13);
    else if (symbol == "Si") return getByAtomicNumber(14);
    else if (symbol == "P")  return getByAtomicNumber(15);
    else if (symbol == "S")  return getByAtomicNumber(16);
    else if (symbol == "Cl") return getByAtomicNumber(17);
    else if (symbol == "Ar") return getByAtomicNumber(18);
    else if (symbol == "K")  return getByAtomicNumber(19);
    else if (symbol == "Ca") return getByAtomicNumber(20);
    else if (symbol == "Sc") return getByAtomicNumber(21);
    else if (symbol == "Ti") return getByAtomicNumber(22);
    else if (symbol == "V")  return getByAtomicNumber(23);
    else if (symbol == "Cr") return getByAtomicNumber(24);
    else if (symbol == "Mn") return getByAtomicNumber(25);
    else if (symbol == "Fe") return getByAtomicNumber(26);
    else if (symbol == "Co") return getByAtomicNumber(27);
    else if (symbol == "Ni") return getByAtomicNumber(28);
    else if (symbol == "Cu") return getByAtomicNumber(29);
    else if (symbol == "Zn") return getByAtomicNumber(30);
    else if (symbol == "Ga") return getByAtomicNumber(31);
    else if (symbol == "Ge") return getByAtomicNumber(32);
    else if (symbol == "As") return getByAtomicNumber(33);
    else if (symbol == "Se") return getByAtomicNumber(34);
    else if (symbol == "Br") return getByAtomicNumber(35);
    else if (symbol == "Kr") return getByAtomicNumber(36);
    else if (symbol == "Rb") return getByAtomicNumber(37);
    else if (symbol == "Sr") return getByAtomicNumber(38);
    else if (symbol == "Y")  return getByAtomicNumber(39);
    else if (symbol == "Zr") return getByAtomicNumber(40);
    else if (symbol == "Nb") return getByAtomicNumber(41);
    else if (symbol == "Mo") return getByAtomicNumber(42);
    else if (symbol == "Tc") return getByAtomicNumber(43);
    else if (symbol == "Ru") return getByAtomicNumber(44);
    else if (symbol == "Rh") return getByAtomicNumber(45);
    else if (symbol == "Pd") return getByAtomicNumber(46);
    else if (symbol == "Ag") return getByAtomicNumber(47);
    else if (symbol == "Cd") return getByAtomicNumber(48);
    else if (symbol == "In") return getByAtomicNumber(49);
    else if (symbol == "Sn") return getByAtomicNumber(50);
    else if (symbol == "Sb") return getByAtomicNumber(51);
    else if (symbol == "Te") return getByAtomicNumber(52);
    else if (symbol == "I")  return getByAtomicNumber(53);
    else if (symbol == "Xe") return getByAtomicNumber(54);
    else if (symbol == "Cs") return getByAtomicNumber(55);
    else if (symbol == "Ba") return getByAtomicNumber(56);
    else if (symbol == "La") return getByAtomicNumber(57);
    else if (symbol == "Ce") return getByAtomicNumber(58);
    else if (symbol == "Pr") return getByAtomicNumber(59);
    else if (symbol == "Nd") return getByAtomicNumber(60);
    else if (symbol == "Pm") return getByAtomicNumber(61);
    else if (symbol == "Sm") return getByAtomicNumber(62);
    else if (symbol == "Eu") return getByAtomicNumber(63);
    else if (symbol == "Gd") return getByAtomicNumber(64);
    else if (symbol == "Tb") return getByAtomicNumber(65);
    else if (symbol == "Dy") return getByAtomicNumber(66);
    else if (symbol == "Ho") return getByAtomicNumber(67);
    else if (symbol == "Er") return getByAtomicNumber(68);
    else if (symbol == "Tm") return getByAtomicNumber(69);
    else if (symbol == "Yb") return getByAtomicNumber(70);
    else if (symbol == "Lu") return getByAtomicNumber(71);
    else if (symbol == "Hf") return getByAtomicNumber(72);
    else if (symbol == "Ta") return getByAtomicNumber(73);
    else if (symbol == "W")  return getByAtomicNumber(74);
    else if (symbol == "Re") return getByAtomicNumber(75);
    else if (symbol == "Os") return getByAtomicNumber(76);
    else if (symbol == "Ir") return getByAtomicNumber(77);
    else if (symbol == "Pt") return getByAtomicNumber(78);
    else if (symbol == "Au") return getByAtomicNumber(79);
    else if (symbol == "Hg") return getByAtomicNumber(80);
    else if (symbol == "Tl") return getByAtomicNumber(81);
    else if (symbol == "Pb") return getByAtomicNumber(82);
    else if (symbol == "Bi") return getByAtomicNumber(83);
    else if (symbol == "Po") return getByAtomicNumber(84);
    else if (symbol == "At") return getByAtomicNumber(85);
    else if (symbol == "Rn") return getByAtomicNumber(86);
    else if (symbol == "Fr") return getByAtomicNumber(87);
    else if (symbol == "Ra") return getByAtomicNumber(88);
    else if (symbol == "Ac") return getByAtomicNumber(89);
    else if (symbol == "Th") return getByAtomicNumber(90);
    else if (symbol == "Pa") return getByAtomicNumber(91);
    else if (symbol == "U")  return getByAtomicNumber(92);
    else if (symbol == "Np") return getByAtomicNumber(93);
    else if (symbol == "Pu") return getByAtomicNumber(94);
    else if (symbol == "Am") return getByAtomicNumber(95);
    else if (symbol == "Cm") return getByAtomicNumber(96);
    else if (symbol == "Bk") return getByAtomicNumber(97);
    else if (symbol == "Cf") return getByAtomicNumber(98);
    else if (symbol == "Es") return getByAtomicNumber(99);
    else if (symbol == "Fm") return getByAtomicNumber(100);
    else if (symbol == "Md") return getByAtomicNumber(101);
    else if (symbol == "No") return getByAtomicNumber(102);
    else if (symbol == "Lr") return getByAtomicNumber(103);
    else if (symbol == "Rf") return getByAtomicNumber(104);
    else if (symbol == "Db") return getByAtomicNumber(105);
    else if (symbol == "Sg") return getByAtomicNumber(106);
    else if (symbol == "Bh") return getByAtomicNumber(107);
    else if (symbol == "Hs") return getByAtomicNumber(108);
    else if (symbol == "Mt") return getByAtomicNumber(109);
    else if (symbol == "Ds") return getByAtomicNumber(110);
    else if (symbol == "Rg") return getByAtomicNumber(111);
    else if (symbol == "Uub" || symbol == "Cn") return getByAtomicNumber(112);
    else if (symbol == "Uut" || symbol == "Nh") return getByAtomicNumber(113);
    else if (symbol == "Uuq" || symbol == "Fl") return getByAtomicNumber(114);
    else if (symbol == "Uup" || symbol == "Mc") return getByAtomicNumber(115);
    else if (symbol == "Uuh" || symbol == "Lv") return getByAtomicNumber(116);

    throw std::out_of_range("Invalid element symbol " + symbol);
}

std::ostream& operator<<(std::ostream& o, const Element& e) {
    return o << e.getName();
}

Element::Element() 
//    : HandleBase(new ElementRep())
    : impl(new Impl())
{}

Element::Element(int atomicNumber, Name name, Symbol symbol, mdunits::Mass typicalMass) 
    : impl(new Impl(atomicNumber, name, symbol, typicalMass))
{}

Element::Element(const Element &src)
    : impl(src.impl)
{}

Element::Element(Element &&src) noexcept
    : impl(std::move(src.impl))
{}

Element::Symbol Element::getSymbol() const {
    return impl->getSymbol();
}
Element::Symbol Element::getName() const {
    return impl->getName();
}
int Element::getAtomicNumber() const {
	return impl->getAtomicNumber();
}

mdunits::Mass Element::getMass() const {
    return impl->getMass();
}

Element & Element::operator=(const Element &src) {
    impl = src.impl;
    return *this;
}

Element & Element::operator=(Element &&src) noexcept {
    impl = std::move(src.impl);
    return *this;
}

//const Element& Element::Hydrogen() {
//    if (elementsByAtomicNumber.find(1) == elementsByAtomicNumber.end())
//        defineCanonicalElement(1, "hydrogen", "H", 1.007947);
//    return elementsByAtomicNumber.find(1)->second;
//}
//const Element& Element::Helium() {
//    if (elementsByAtomicNumber.find(2) == elementsByAtomicNumber.end())
//        defineCanonicalElement(2, "helium", "He", 4.003);
//    return elementsByAtomicNumber.find(2)->second;
//}
//const Element& Element::Carbon() {
//    if (elementsByAtomicNumber.find(6) == elementsByAtomicNumber.end())
//        defineCanonicalElement(6, "carbon", "C", 12.01078);
//    return elementsByAtomicNumber.find(6)->second;
//}
//const Element& Element::Nitrogen() {
//    if (elementsByAtomicNumber.find(7) == elementsByAtomicNumber.end())
//        defineCanonicalElement(7, "nitrogen", "N", 14.00672);
//    return elementsByAtomicNumber.find(7)->second;
//}
//const Element& Element::Oxygen() {
//    if (elementsByAtomicNumber.find(8) == elementsByAtomicNumber.end())
//        defineCanonicalElement(8, "oxygen", "O", 15.99943);
//    return elementsByAtomicNumber.find(8)->second;
//}
//const Element& Element::Neon() {
//    if (elementsByAtomicNumber.find(10) == elementsByAtomicNumber.end())
//        defineCanonicalElement(10, "neon", "Ne", 20.180);
//    return elementsByAtomicNumber.find(10)->second;
//}
//const Element& Element::Phosphorus() {
//    if (elementsByAtomicNumber.find(15) == elementsByAtomicNumber.end())
//        defineCanonicalElement(15, "phosphorus", "P", 30.9737622);
//    return elementsByAtomicNumber.find(15)->second;
//}
//const Element& Element::Sulfur() {
//    if (elementsByAtomicNumber.find(16) == elementsByAtomicNumber.end())
//        defineCanonicalElement(16, "sulfur", "S", 32.0655);
//    return elementsByAtomicNumber.find(16)->second;
//}
//const Element& Element::Argon() {
//    if (elementsByAtomicNumber.find(18) == elementsByAtomicNumber.end())
//        defineCanonicalElement(18, "argon", "Ar", 39.9481);
//    return elementsByAtomicNumber.find(18)->second;
//}
//const Element& Element::Selenium() {
//    if (elementsByAtomicNumber.find(34) == elementsByAtomicNumber.end())
//        defineCanonicalElement(34, "selenium", "Se", 78.963);
//    return elementsByAtomicNumber.find(34)->second;
//}
//const Element& Element::Xenon() {
//    if (elementsByAtomicNumber.find(54) == elementsByAtomicNumber.end())
//        defineCanonicalElement(54, "xenon", "Xe", 131.290);
//    return elementsByAtomicNumber.find(54)->second;
//}
//const Element& Element::Radon() {
//    if (elementsByAtomicNumber.find(86) == elementsByAtomicNumber.end())
//        defineCanonicalElement(86, "radon", "Rn", 222.018);
//    return elementsByAtomicNumber.find(86)->second;
//}

} // namespace SimTK
