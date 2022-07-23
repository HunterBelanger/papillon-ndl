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
#include <string>

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

    if (A_ > 300 && A_ < 600) {
      A_ -= 300;
      level_ = 1;
    } else if (A_ > 600 && A_ < 900) {
      A_ -= 600;
      level_ = 2;
    } else if (A_ > 900) {
      std::string mssg = "ZAID with A > 900 indicates an isomer level > 2. ";
      mssg += "Cannot create a Nuclde with isomer level greater than 2.";
      throw PNDLException(mssg);
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
  ZAID zaid() const { return ZAID(this->Z(), this->A() + 300 * level_); }

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
    return h(nuc.Z() * 1000 + nuc.A() + nuc.level() * 300);
  }
};

#endif
