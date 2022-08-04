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
#ifndef PAPILLON_NDL_ISOTOPE_H
#define PAPILLON_NDL_ISOTOPE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/element.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <ostream>
#include <regex>
#include <string>

namespace pndl {

/**
 * @brief Class which identifies an isotope. The atomic mass number must be at
 *        least equal to the atomic number, and can be no larger than 300.
 */
class Isotope {
 public:
  /**
   * @param element Element of the Isotope.
   * @param A Atomic mass of the Isotope.
   */
  Isotope(const Element& element, uint32_t A) : element_(element), A_(A) {
    if (A < element_.Z()) {
      std::string mssg = "Cannot create isotope " + element_.name() + "-" +
                         std::to_string(this->A()) +
                         ". Isotopes must satisfy A >= Z. "
                         "Was provided with A = " +
                         std::to_string(this->A()) +
                         ", Z = " + std::to_string(this->Z()) + ".";
      throw PNDLException(mssg);
    }

    if (A >= 300) {
      std::string mssg = "Cannot create isotope " + element_.name() + "-" +
                         std::to_string(this->A()) +
                         ". Isotopes must satisfy A < 300."
                         " Was provided with A = " +
                         std::to_string(this->A()) + ".";
      throw PNDLException(mssg);
    }
  }

  /**
   * @param Z Atomic number of the Isotope.
   * @param A Atomic mass of the Isotope.
   */
  Isotope(uint8_t Z, uint32_t A) try : element_(Z), A_(A) {
    if (A_ < element_.Z()) {
      std::string mssg = "Cannot create isotope " + element_.name() + "-" +
                         std::to_string(this->A()) +
                         ". Isotopes must satisfy A >= Z. "
                         "Was provided with A = " +
                         std::to_string(this->A()) +
                         ", Z = " + std::to_string(this->Z()) + ".";
      throw PNDLException(mssg);
    }

    if (A_ >= 300) {
      std::string mssg = "Cannot create isotope " + element_.name() + "-" +
                         std::to_string(this->A()) +
                         ". Isotopes must satisfy A < 300."
                         " Was provided with A = " +
                         std::to_string(this->A()) + ".";
      throw PNDLException(mssg);
    }
  } catch (PNDLException& err) {
    std::string mssg = "Could not construct Element associated with Isotope.";
    err.add_to_exception(mssg);
    throw err;
  }

  /**
   * @param zaid ZAID of the Isotope.
   */
  Isotope(const ZAID& zaid) try : element_(zaid), A_(zaid.A()) {
    if (A_ < element_.Z()) {
      std::string mssg = "Cannot create isotope " + element_.name() + "-" +
                         std::to_string(this->A()) +
                         ". Isotopes must satisfy A >= Z. "
                         "Was provided with A = " +
                         std::to_string(this->A()) +
                         ", Z = " + std::to_string(this->Z()) + ".";
      throw PNDLException(mssg);
    }

    if (A_ >= 300) {
      std::string mssg = "Cannot create isotope " + element_.name() + "-" +
                         std::to_string(this->A()) +
                         ". Isotopes must satisfy A < 300."
                         " Was provided with A = " +
                         std::to_string(this->A()) + ".";
      throw PNDLException(mssg);
    }
  } catch (PNDLException& err) {
    std::string mssg = "Could not construct Element associated with ZAID.";
    err.add_to_exception(mssg);
    throw err;
  }

  /**
   * @param symbol String containing the symbol for the isotope. The symbol
   *               must be in SSAAA format, such as Al27, or U235.
   */
  Isotope(const std::string& symbol) : element_(1), A_(0) {
    const std::regex is_isotope_regex(
        "(^\\s+)?([A-Z][a-z]{0,1}[0-9]{1,3})(\\s+)?");
    if (std::regex_match(symbol, is_isotope_regex) == false) {
      std::string mssg = "The symbol \"" + symbol + "\" is not a valid ";
      mssg += "isotope symbol.";
      throw PNDLException(mssg);
    }

    // We are a valid formatted isotope. First get the element.
    const std::regex element_regex("([A-Z][a-z]{0,1})");
    std::smatch match;
    std::regex_search(symbol, match, element_regex);
    std::string element_symbol(match[0].first, match[0].second);
    try {
      element_ = Element(element_symbol);
    } catch (PNDLException& err) {
      std::string mssg = "Could not create isotope with element symbol \"";
      mssg += element_symbol + "\".";
      err.add_to_exception(mssg);
      throw err;
    }

    // Get the atomic mass number
    const std::regex atomic_mass_regex("([0-9]{1,3})");
    std::regex_search(symbol, match, atomic_mass_regex);
    std::string atomic_mass_str(match[0].first, match[0].second);
    A_ = std::stoul(atomic_mass_str);

    if (A_ < element_.Z()) {
      std::string mssg = "Cannot create isotope " + element_.name() + "-" +
                         std::to_string(this->A()) +
                         ". Isotopes must satisfy A >= Z. "
                         "Was provided with A = " +
                         std::to_string(this->A()) +
                         ", Z = " + std::to_string(this->Z()) + ".";
      throw PNDLException(mssg);
    }

    if (A_ >= 300) {
      std::string mssg = "Cannot create isotope " + element_.name() + "-" +
                         std::to_string(this->A()) +
                         ". Isotopes must satisfy A < 300."
                         " Was provided with A = " +
                         std::to_string(this->A()) + ".";
      throw PNDLException(mssg);
    }
  }

  /**
   * @brief Returns atomic number of isotope.
   */
  uint8_t Z() const { return element_.Z(); }

  /**
   * @brief Returns atomic number of isotope.
   */
  uint8_t atomic_number() const { return Z(); }

  /**
   * @brief Returns atomic mass of isotope.
   */
  uint32_t A() const { return A_; }

  /**
   * @brief Returns atomic mass of isotope.
   */
  uint32_t atomic_mass() const { return A(); }

  /**
   * @brief Returns ZAID of isotope.
   */
  ZAID zaid() const { return ZAID(element_.Z(), A_); }

  /**
   * @brief Returns the symbol of the isotope.
   */
  std::string symbol() const { return element_.symbol() + std::to_string(A_); }

  /**
   * @brief Returns the Element symbol of the isotope.
   */
  const std::string& element_symbol() const { return element_.symbol(); }

  /**
   * @brief Returns the Element name of the isotope.
   */
  const std::string& element_name() const { return element_.name(); }

  /**
   * @brief Returns true if two isotopes are the same, and false if not.
   */
  bool operator==(const Isotope& rhs) const {
    if (Z() == rhs.Z() && A() == rhs.A()) return true;
    return false;
  }

  /**
   * @brief Returns true if one Isotope's atomic number is less than the
   *        other's. If the atomic numbers are equal and the Isotope's atomic
   *        mass is less than the other's, true is also returned. Otherwise,
   *        flase is returned.
   */
  bool operator<(const Isotope& rhs) const {
    if (Z() < rhs.Z())
      return true;
    else if (Z() > rhs.Z())
      return false;

    if (A() < rhs.A()) return true;

    return false;
  }

 private:
  Element element_;
  uint32_t A_;
};

inline std::ostream& operator<<(std::ostream& strm, const Isotope& istp) {
  strm << istp.element_symbol() << istp.A();
  return strm;
}

}  // namespace pndl

template <>
struct std::hash<pndl::Isotope> {
  std::size_t operator()(const pndl::Isotope& iso) const noexcept {
    std::hash<uint32_t> h;
    return h(iso.Z() * 1000 + iso.A());
  }
};

#endif
