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
#ifndef PAPILLON_NDL_ELEMENT_H
#define PAPILLON_NDL_ELEMENT_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/zaid.hpp>
#include <array>
#include <cstdint>
#include <functional>
#include <ostream>
#include <string>

namespace pndl {

/**
 * @brief Class which identifies an element.
 */
class Element {
 public:
  /**
   * @param Z Atomic number of the element. Must be in the interval [1,118].
   */
  Element(uint8_t Z) : Z_(Z) {
    if (Z_ == 0 || Z_ > N_ELEM) {
      std::string mssg =
          "Elements must have an atomic number "
          "in interval [1," +
          std::to_string(N_ELEM) + "].";
      throw PNDLException(mssg);
    }
  }

  /**
   * @param zaid ZAID identifier. Requires zaid.Z() be in the interval [1,118].
   */
  Element(const ZAID& zaid) : Z_(zaid.Z()) {
    if (Z_ == 0 || Z_ > N_ELEM) {
      std::string mssg =
          "Elements must have an atomic number "
          "in interval [1," +
          std::to_string(N_ELEM) + "].";
      throw PNDLException(mssg);
    }
  }

  /**
   * @brief Returns the atomic number of the element.
   */
  uint8_t Z() const { return Z_; }

  /**
   * @brief Returns the atomic number of the element.
   */
  uint8_t atomic_number() const { return Z_; }

  /**
   * @brief Returns the symbol of the element.
   */
  const std::string& symbol() const {
    return elements_table[static_cast<std::size_t>(Z_ - 1)].symbol;
  }

  /**
   * @brief Returns the name of the element.
   */
  const std::string& name() const {
    return elements_table[static_cast<std::size_t>(Z_ - 1)].name;
  }

  /**
   * @brief Returns the ZAID which represents the natural element.
   */
  ZAID zaid() const { return ZAID(Z_, 0); }

  /**
   * @brief Returns true if two elements are the same, and false if not.
   */
  bool operator==(const Element& rhs) const { return this->Z_ == rhs.Z_; }

  /**
   * @brief Returns true if the Element's atomic number is less than the
   *        other's, and false if otherwise.
   */
  bool operator<(const Element& rhs) const { return this->Z_ < rhs.Z_; }

  /**
   * @brief Finds an element from a symbol.
   * @param symbol String which holds the element symbol for which to search.
   */
  static Element from_symbol(const std::string& symbol);

  /**
   * @brief Finds an element from a name.
   * @param symbol String which holds the element name for which to search.
   */
  static Element from_name(const std::string& name);

 private:
  struct Info {
    std::string name;
    std::string symbol;
  };
  static constexpr uint8_t N_ELEM{118};

  uint8_t Z_;

  static std::array<Info, N_ELEM> elements_table;
};

inline std::ostream& operator<<(std::ostream& strm, const Element& elem) {
  strm << elem.symbol();
  return strm;
}

}  // namespace pndl

template <>
struct std::hash<pndl::Element> {
  std::size_t operator()(const pndl::Element& elem) const noexcept {
    std::hash<uint8_t> h;
    return h(elem.Z());
  }
};

#endif
