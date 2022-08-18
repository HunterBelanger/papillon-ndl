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
#ifndef PAPILLON_NDL_NUCLIDE_H
#define PAPILLON_NDL_NUCLIDE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/isotope.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/zaid.hpp>
#include <cstdint>
#include <functional>
#include <ostream>
#include <regex>
#include <string>
#include <sstream>

namespace pndl {

/**
 * @brief Class which identifies a nuclide. The isomer level may be no
 *        greater than 2.
 */
class Nuclide {
 public:
  /**
   * @param isotope Isotope of the nuclide.
   * @param level Isomer level of the nuclide.
   */
  Nuclide(const Isotope& isotope, uint8_t level = 0)
      : isotope_(isotope), level_(level) {
    if (level > 2) {
      std::string mssg =
          "Cannot create Nuclide with iosmer lever greater than 2. "
          "Was provided with level = " +
          std::to_string(level_) + ".";
      throw PNDLException(mssg);
    }
  }

  /**
   * @param Z Atomic number of nuclide.
   * @param A Atomic mass of nuclide.
   * @param level Isomer level of nuclide.
   */
  Nuclide(uint8_t Z, uint32_t A, uint8_t level = 0) try
      : isotope_(Z, A), level_(level) {
    if (level > 2) {
      std::string mssg =
          "Cannot create Nuclide with iosmer lever greater than 2. "
          "Was provided with level = " +
          std::to_string(level_) + ".";
      throw PNDLException(mssg);
    }
  } catch (PNDLException& err) {
    std::string mssg = "Could not create isotope.";
    err.add_to_exception(mssg);
    throw err;
  }

  /**
   * @param zaid ZAID of nuclide.
   */
  Nuclide(const ZAID& zaid) : isotope_(1, 1), level_(0) {
    uint8_t Z_ = zaid.Z();
    uint32_t A_ = zaid.A();
    level_ = 0;
    Element el(1);
    try {
      el = Element(Z_);
    } catch (PNDLException& err) {
      std::string mssg = "Could not create Nuclide for Z = " + std::to_string(Z_) + ".";
      err.add_to_exception(mssg);
      throw err;
    }

    if (A_ > 300) {
      level_++;
      A_ -= 400;

      while (A_ > el.largest_isotope() && A_ > 100) {
        A_ -= 100;
        level_++;
      }

      if (A_ > el.largest_isotope()) {
        std::stringstream mssg;
        mssg << "ZAID = " << zaid << " indicates Z = " << +Z_ << ", A = " << A_;
        mssg << ", m = " << +level_ << ". The largest possible atomic mass for ";
        mssg << el.symbol() << " is " << el.largest_isotope() << ".";
        throw PNDLException(mssg.str());
      }
    }

    try {
      isotope_ = Isotope(Z_, A_);
    } catch (PNDLException& err) {
      std::string mssg = "Could not create isotope.";
      err.add_to_exception(mssg);
      throw err;
    }
  }

  /**
   * @param symbol String containing the symbol for the nuclide. The symbol
   *               must be in SSAAA format. If the nuclide is an isomer, the
   *               isomer level can be added as SSAAAmL. The level can be
   *               0, 1, or 2.
   */
  Nuclide(const std::string& symbol) : isotope_(1, 1), level_(0) {
    const std::regex is_nuclide_regex(
        "(^\\s+)?([A-Z][a-z]{0,1}[0-9]{1,3})([m][0-2])?(\\s+)?");

    if (std::regex_match(symbol, is_nuclide_regex) == false) {
      std::string mssg = "The symbol \"" + symbol + "\" is not a valid ";
      mssg += "Nuclide symbol.";
      throw PNDLException(mssg);
    }

    const std::regex isotope_regex("([A-Z][a-z]{0,1}[0-9]{1,3})");
    std::smatch match;
    std::regex_search(symbol, match, isotope_regex);
    std::string isotope_symbol(match[0].first, match[0].second);
    try {
      isotope_ = Isotope(isotope_symbol);
    } catch (PNDLException& err) {
      std::string mssg = "Could not create nuclide with isotope symbol \"";
      mssg += isotope_symbol + "\".";
      err.add_to_exception(mssg);
      throw err;
    }

    const std::regex isomer_regex("([m][0-2])");
    std::regex_search(symbol, match, isomer_regex);
    std::string isomer_str(match[0].first, match[0].second);
    level_ = 0;
    if (isomer_str.size() > 0) {
      isomer_str.erase(isomer_str.begin());
      level_ = static_cast<uint8_t>(std::stoul(isomer_str));
    }
  }

  /**
   * @brief Returns atomic number of isotope.
   */
  uint8_t Z() const { return isotope_.Z(); }

  /**
   * @brief Returns atomic number of isotope.
   */
  uint8_t atomic_number() const { return this->Z(); }

  /**
   * @brief Returns atomic mass of isotope.
   */
  uint32_t A() const { return isotope_.A(); }

  /**
   * @brief Returns atomic mass of isotope.
   */
  uint32_t atomic_mass() const { return this->A(); }

  /**
   * @brief Returns the isomer level of the nuclide.
   */
  uint8_t level() const { return level_; }

  /**
   * @brief Returns the ZAID for the nuclide.
   */
  ZAID zaid() const {
    uint32_t A_za = this->A();
    if (this->level() > 0) {
      A_za += 300;
      A_za += static_cast<uint32_t>(this->level()) * 100;
    }
    return ZAID(this->Z(), A_za);
  }

  /**
   * @brief Returns the symbol of the nuclide.
   */
  std::string symbol() const {
    std::string symbl = isotope_.symbol();
    if (level_ > 0) {
      symbl += 'm' + std::to_string(level_);
    }
    return symbl;
  }

  /**
   * @brief Returns the symbol of the isotope.
   */
  std::string isotope_symbol() const { return isotope_.symbol(); }

  /**
   * @brief Returns the element symbol of the nuclide.
   */
  const std::string& element_symbol() const {
    return isotope_.element_symbol();
  }

  /**
   * @brief Returns the element name of the nuclide.
   */
  const std::string& element_name() const { return isotope_.element_name(); }

  /**
   * @brief Returns true if two Nuclide are the same, and false if not.
   */
  bool operator==(const Nuclide& rhs) const {
    if (Z() == rhs.Z() && A() == rhs.A() && level() == rhs.level()) return true;
    return false;
  }

  /**
   * @brief Returns true if one Nuclides's atomic number is less than the
   *        other's. If the atomic numbers are equal and the Nuclides's atomic
   *        mass is less than the other's, true is also returned. If the
   *        atomic masses are equal, and the isomer level is less than the
   *        other's, true is returned. Otherwise, flase is returned.
   */
  bool operator<(const Nuclide& rhs) const {
    if (Z() < rhs.Z())
      return true;
    else if (Z() > rhs.Z())
      return false;

    if (A() < rhs.A())
      return true;
    else if (A() > rhs.A())
      return false;

    if (level() < rhs.level()) return true;

    return false;
  }

 private:
  Isotope isotope_;
  uint8_t level_;
};

inline std::ostream& operator<<(std::ostream& strm, const Nuclide& nuc) {
  strm << nuc.symbol();
  return strm;
}

}  // namespace pndl

template <>
struct std::hash<pndl::Nuclide> {
  std::size_t operator()(const pndl::Nuclide& nuc) const noexcept {
    std::hash<uint32_t> h;
    uint32_t A = nuc.A();
    if (nuc.level() > 0) {
      A += 300;
      A += static_cast<uint32_t>(nuc.level()) * 100;
    }
    return h(nuc.Z() * 1000 + A);
  }
};

#endif
