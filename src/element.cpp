/*
 * Papillon Nuclear Data Library
 * Copyright 2021, Hunter Belanger
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

namespace pndl {

Element Element::from_symbol(const std::string& symbol) {
  uint8_t z = 0;
  bool found = false;
  for (z = 0; z < elements_table.size(); z++) {
    if (elements_table[z].symbol == symbol) {
      found = true;
      break; 
    } 
  }
  
  // Make sure we found the element
  if (found == false) {
    std::string mssg = "Could not find an element with symbol \"" + symbol + "\"."; 
    throw PNDLException(mssg);
  }
  
  // Advance z by 1 for the correct atomic number.
  z++;

  return Element(z);
}

Element Element::from_name(const std::string& name) {
  uint8_t z = 0;
  bool found = false;
  for (z = 0; z < elements_table.size(); z++) {
    if (elements_table[z].name == name) {
      found = true;
      break; 
    } 
  }
  
  // Make sure we found the element
  if (found == false) {
    std::string mssg = "Could not find an element with name \"" + name + "\"."; 
    throw PNDLException(mssg);
  }
  
  // Advance z by 1 for the correct atomic number.
  z++;

  return Element(z);
}

std::array<Element::Info, Element::N_ELEM> Element::elements_table {
  Element::Info {"Hydrogen", "H"},      Element::Info {"Helium", "He"},
  Element::Info {"Lithium", "Li"},      Element::Info {"Beryllium", "Be"},
  Element::Info {"Boron", "B"},         Element::Info {"Carbon", "C"},
  Element::Info {"Nitrogen", "N"},      Element::Info {"Oxygen", "O"},
  Element::Info {"Flourine", "F"},      Element::Info {"Neon", "Ne"},
  Element::Info {"Sodium", "Na"},       Element::Info {"Magnesium", "Mg"},
  Element::Info {"Aluminum", "Al"},     Element::Info {"Silicon", "Si"},
  Element::Info {"Phosphorus", "P"},    Element::Info {"Sulfur", "S"},
  Element::Info {"Chlorine", "Cl"},     Element::Info {"Argon", "Ar"},
  Element::Info {"Potassium", "K"},     Element::Info {"Calcium", "Ca"},
  Element::Info {"Scandium", "Sc"},     Element::Info {"Titanium", "Ti"},
  Element::Info {"Vanadium", "V"},      Element::Info {"Chromium", "Cr"},
  Element::Info {"Manganese", "Mn"},    Element::Info {"Iron", "Fe"},
  Element::Info {"Cobalt", "Co"},       Element::Info {"Nickel", "Ni"},
  Element::Info {"Copper", "Cu"},       Element::Info {"Zinc", "Zn"},
  Element::Info {"Gallium", "Ga"},      Element::Info {"Germanium", "Ge"},
  Element::Info {"Arsenic", "As"},      Element::Info {"Selenium", "Se"},
  Element::Info {"Bromine", "Br"},      Element::Info {"Krypton", "Kr"},
  Element::Info {"Rubidium", "Rb"},     Element::Info {"Strontium", "Sr"},
  Element::Info {"Yttrium", "Y"},       Element::Info {"Zirconium", "Zr"},
  Element::Info {"Niobium", "Nb"},      Element::Info {"Molbdenum", "Mo"},
  Element::Info {"Technetium", "Tc"},   Element::Info {"Ruthenium", "Ru"},
  Element::Info {"Rhodium", "Rh"},      Element::Info {"Palladium", "Pd"},
  Element::Info {"Silver", "Ag"},       Element::Info {"Cadmium", "Cd"},
  Element::Info {"Indium", "In"},       Element::Info {"Tin", "Sn"},
  Element::Info {"Antimony", "Sb"},     Element::Info {"Tellurium", "Te"},
  Element::Info {"Iodine", "I"},        Element::Info {"Xenon", "Xe"},
  Element::Info {"Cesium", "Cs"},       Element::Info {"Barium", "Ba"},
  Element::Info {"Lanthanum", "La"},    Element::Info {"Cerium", "Ce"},
  Element::Info {"Praseodymium", "Pr"}, Element::Info {"Neodymium", "Nd"},
  Element::Info {"Promethium", "Pm"},   Element::Info {"Samarium", "Sm"},
  Element::Info {"Europium", "Eu"},     Element::Info {"Gadolinium", "Gd"},
  Element::Info {"Terbium", "Tb"},      Element::Info {"Dysprosium", "Dy"},
  Element::Info {"Holmium", "Ho"},      Element::Info {"Erbium", "Er"},
  Element::Info {"Thulium", "Tm"},      Element::Info {"Ytterbium", "Yb"},
  Element::Info {"Lutetium", "Lu"},     Element::Info {"Hafnium", "Hf"},
  Element::Info {"Tantalum", "Ta"},     Element::Info {"Tungsten", "W"},
  Element::Info {"Rhenium", "Re"},      Element::Info {"Osmium", "Os"},
  Element::Info {"Iridium", "Ir"},      Element::Info {"Platinum", "Pt"},
  Element::Info {"Gold", "Au"},         Element::Info {"Mercury", "Hg"},
  Element::Info {"Thallium", "Tl"},     Element::Info {"Lead", "Pb"},
  Element::Info {"Bismuth", "Bi"},      Element::Info {"Polonium", "Po"},
  Element::Info {"Astatine", "At"},     Element::Info {"Radon", "Rn"},
  Element::Info {"Francium", "Fr"},     Element::Info {"Radium", "Ra"},
  Element::Info {"Actinium", "Ac"},     Element::Info {"Thorium", "Th"},
  Element::Info {"Protactinium", "Pa"}, Element::Info {"Uranium", "U"},
  Element::Info {"Neptunium", "Np"},    Element::Info {"Plutonium", "Pu"},
  Element::Info {"Americium", "Am"},    Element::Info {"Curium", "Cm"},
  Element::Info {"Berkelium", "Bk"},    Element::Info {"Californium", "Cf"},
  Element::Info {"Einsteinium", "Es"},  Element::Info {"Fermium", "Fm"},
  Element::Info {"Mendelevium", "Md"},  Element::Info {"Nobelium", "No"},
  Element::Info {"Lawrencium", "Lr"},   Element::Info {"Rutherfordium", "Rf"},
  Element::Info {"Dubnium", "Db"},      Element::Info {"Seaborgium", "Sg"},
  Element::Info {"Bohrium", "Bh"},      Element::Info {"Hassium", "Hs"},
  Element::Info {"Meitnerium", "Mt"},   Element::Info {"Darmstadtium", "Ds"},
  Element::Info {"Roentgenium", "Rg"},  Element::Info {"Copernicium", "Cn"},
  Element::Info {"Nihonium", "Nh"},     Element::Info {"Flerovium", "Fl"},
  Element::Info {"Moscovium", "Mc"},    Element::Info {"Livermorium", "Lv"},
  Element::Info {"Tennessine", "Ts"},   Element::Info {"Oganesson", "Og"}
};

}
