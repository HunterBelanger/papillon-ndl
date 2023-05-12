/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2023, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Papillon Nuclear Data Library (PapillonNDL).
 *
 * PapillonNDL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PapillonNDL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PapillonNDL. If not, see <https://www.gnu.org/licenses/>.
 *
 * */

#include <PapillonNDL/element.hpp>
#include <regex>

namespace pndl {

Element::Element(const std::string& name_or_symbol) : Z_(0) {
  bool found = false;

  // First check and see if we are a valid symbol.
  const std::regex is_element_regex("(^\\s+)?([A-Z][a-z]{0,1})(\\s+)?");
  if (std::regex_match(name_or_symbol, is_element_regex)) {
    const std::regex element_regex("([A-Z][a-z]{0,1})");
    std::smatch match;
    std::regex_search(name_or_symbol, match, element_regex);
    std::string element_symbol(match[0].first, match[0].second);

    for (Z_ = 0; Z_ < elements_table.size(); Z_++) {
      if (elements_table[Z_].symbol == element_symbol) {
        found = true;
        break;
      }
    }
  }

  // If we still haven't found it, maybe it's a name
  const std::regex is_name_regex("(^\\s+)?\\b([A-Z][a-z]+)\\b(\\s+)?");
  if (found == false && std::regex_match(name_or_symbol, is_name_regex)) {
    Z_ = 0;
    const std::regex name_regex("([A-Z][a-z]+)");
    std::smatch match;
    std::regex_search(name_or_symbol, match, name_regex);
    std::string element_name(match[0].first, match[0].second);

    for (Z_ = 0; Z_ < N_ELEM; Z_++) {
      if (elements_table[Z_].name == element_name) {
        found = true;
        break;
      }
    }
  }

  if (found) {
    // Advance z by 1 for the correct atomic number.
    Z_++;
  } else {
    // We couldn't find the element...
    std::string mssg = "Could not find an element symbol or name matching \"";
    mssg += name_or_symbol + "\".";
    throw PNDLException(mssg);
  }
}

std::array<Element::Info, Element::N_ELEM> Element::elements_table{
    Element::Info{"Hydrogen", "H", 7},
    Element::Info{"Helium", "He", 10},
    Element::Info{"Lithium", "Li", 13},
    Element::Info{"Beryllium", "Be", 16},
    Element::Info{"Boron", "B", 21},
    Element::Info{"Carbon", "C", 23},
    Element::Info{"Nitrogen", "N", 25},
    Element::Info{"Oxygen", "O", 28},
    Element::Info{"Flourine", "F", 31},
    Element::Info{"Neon", "Ne", 34},
    Element::Info{"Sodium", "Na", 37},
    Element::Info{"Magnesium", "Mg", 40},
    Element::Info{"Aluminum", "Al", 43},
    Element::Info{"Silicon", "Si", 45},
    Element::Info{"Phosphorus", "P", 47},
    Element::Info{"Sulfur", "S", 49},
    Element::Info{"Chlorine", "Cl", 51},
    Element::Info{"Argon", "Ar", 53},
    Element::Info{"Potassium", "K", 56},
    Element::Info{"Calcium", "Ca", 58},
    Element::Info{"Scandium", "Sc", 61},
    Element::Info{"Titanium", "Ti", 63},
    Element::Info{"Vanadium", "V", 66},
    Element::Info{"Chromium", "Cr", 68},
    Element::Info{"Manganese", "Mn", 71},
    Element::Info{"Iron", "Fe", 74},
    Element::Info{"Cobalt", "Co", 76},
    Element::Info{"Nickel", "Ni", 79},
    Element::Info{"Copper", "Cu", 82},
    Element::Info{"Zinc", "Zn", 85},
    Element::Info{"Gallium", "Ga", 87},
    Element::Info{"Germanium", "Ge", 90},
    Element::Info{"Arsenic", "As", 92},
    Element::Info{"Selenium", "Se", 95},
    Element::Info{"Bromine", "Br", 98},
    Element::Info{"Krypton", "Kr", 101},
    Element::Info{"Rubidium", "Rb", 103},
    Element::Info{"Strontium", "Sr", 107},
    Element::Info{"Yttrium", "Y", 109},
    Element::Info{"Zirconium", "Zr", 112},
    Element::Info{"Niobium", "Nb", 115},
    Element::Info{"Molbdenum", "Mo", 117},
    Element::Info{"Technetium", "Tc", 120},
    Element::Info{"Ruthenium", "Ru", 124},
    Element::Info{"Rhodium", "Rh", 126},
    Element::Info{"Palladium", "Pd", 128},
    Element::Info{"Silver", "Ag", 130},
    Element::Info{"Cadmium", "Cd", 133},
    Element::Info{"Indium", "In", 135},
    Element::Info{"Tin", "Sn", 138},
    Element::Info{"Antimony", "Sb", 140},
    Element::Info{"Tellurium", "Te", 143},
    Element::Info{"Iodine", "I", 145},
    Element::Info{"Xenon", "Xe", 148},
    Element::Info{"Cesium", "Cs", 151},
    Element::Info{"Barium", "Ba", 153},
    Element::Info{"Lanthanum", "La", 155},
    Element::Info{"Cerium", "Ce", 157},
    Element::Info{"Praseodymium", "Pr", 159},
    Element::Info{"Neodymium", "Nd", 161},
    Element::Info{"Promethium", "Pm", 163},
    Element::Info{"Samarium", "Sm", 165},
    Element::Info{"Europium", "Eu", 167},
    Element::Info{"Gadolinium", "Gd", 169},
    Element::Info{"Terbium", "Tb", 171},
    Element::Info{"Dysprosium", "Dy", 173},
    Element::Info{"Holmium", "Ho", 175},
    Element::Info{"Erbium", "Er", 177},
    Element::Info{"Thulium", "Tm", 179},
    Element::Info{"Ytterbium", "Yb", 181},
    Element::Info{"Lutetium", "Lu", 185},
    Element::Info{"Hafnium", "Hf", 189},
    Element::Info{"Tantalum", "Ta", 192},
    Element::Info{"Tungsten", "W", 194},
    Element::Info{"Rhenium", "Re", 198},
    Element::Info{"Osmium", "Os", 202},
    Element::Info{"Iridium", "Ir", 204},
    Element::Info{"Platinum", "Pt", 206},
    Element::Info{"Gold", "Au", 210},
    Element::Info{"Mercury", "Hg", 216},
    Element::Info{"Thallium", "Tl", 218},
    Element::Info{"Lead", "Pb", 220},
    Element::Info{"Bismuth", "Bi", 224},
    Element::Info{"Polonium", "Po", 227},
    Element::Info{"Astatine", "At", 229},
    Element::Info{"Radon", "Rn", 231},
    Element::Info{"Francium", "Fr", 233},
    Element::Info{"Radium", "Ra", 235},
    Element::Info{"Actinium", "Ac", 237},
    Element::Info{"Thorium", "Th", 239},
    Element::Info{"Protactinium", "Pa", 241},
    Element::Info{"Uranium", "U", 243},
    Element::Info{"Neptunium", "Np", 245},
    Element::Info{"Plutonium", "Pu", 247},
    Element::Info{"Americium", "Am", 249},
    Element::Info{"Curium", "Cm", 252},
    Element::Info{"Berkelium", "Bk", 254},
    Element::Info{"Californium", "Cf", 256},
    Element::Info{"Einsteinium", "Es", 258},
    Element::Info{"Fermium", "Fm", 260},
    Element::Info{"Mendelevium", "Md", 262},
    Element::Info{"Nobelium", "No", 264},
    Element::Info{"Lawrencium", "Lr", 266},
    Element::Info{"Rutherfordium", "Rf", 268},
    Element::Info{"Dubnium", "Db", 270},
    Element::Info{"Seaborgium", "Sg", 273},
    Element::Info{"Bohrium", "Bh", 275},
    Element::Info{"Hassium", "Hs", 277},
    Element::Info{"Meitnerium", "Mt", 279},
    Element::Info{"Darmstadtium", "Ds", 281},
    Element::Info{"Roentgenium", "Rg", 283},
    Element::Info{"Copernicium", "Cn", 285},
    Element::Info{"Nihonium", "Nh", 287},
    Element::Info{"Flerovium", "Fl", 289},
    Element::Info{"Moscovium", "Mc", 291},
    Element::Info{"Livermorium", "Lv", 293},
    Element::Info{"Tennessine", "Ts", 294},
    Element::Info{"Oganesson", "Og", 295}};

}  // namespace pndl
