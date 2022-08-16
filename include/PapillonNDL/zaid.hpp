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
#ifndef PAPILLON_NDL_ZAID_H
#define PAPILLON_NDL_ZAID_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <cstdint>
#include <functional>
#include <ostream>

namespace pndl {

/**
 * @brief Class which represents a ZAID identifier.
 */
class ZAID {
 public:
  /**
   * @param Z Atomic number for ZAID.
   * @param A Atomic mass for ZAID.
   */
  ZAID(uint8_t Z, uint32_t A) : Z_(Z), A_(A) {}

  /**
   * @brief Returns the atomic number of ZAID.
   */
  uint8_t Z() const { return Z_; }

  /**
   * @brief Returns the atomic mass of ZAID.
   */
  uint32_t A() const { return A_; }

  /**
   * @brief Returns the ZAID as an unsigned integer.
   */
  uint32_t zaid() const { return 1000 * Z_ + A_; };

  /**
   * @brief Returns true if two ZAIDs are equal, and false if not.
   */
  bool operator==(const ZAID& rhs) const {
    return (this->Z_ == rhs.Z_) && (this->A_ == rhs.A_);
  }

  /**
   * @brief Returns true if one ZAID's atomic number is less than the
   *        other's. If the atomic numbers are equal and the ZAID's atomic
   *        mass is less than the other's, true is also returned. Otherwise,
   *        flase is returned.
   */
  bool operator<(const ZAID& rhs) const {
    if (Z() < rhs.Z())
      return true;
    else if (Z() > rhs.Z())
      return false;

    if (A() < rhs.A()) return true;

    return false;
  }

 private:
  uint8_t Z_;
  uint32_t A_;
};

inline std::ostream& operator<<(std::ostream& strm, const ZAID& zaid) {
  strm << +zaid.Z() << zaid.A();
  return strm;
}

}  // namespace pndl

template <>
struct std::hash<pndl::ZAID> {
  std::size_t operator()(const pndl::ZAID& zaid) const noexcept {
    std::hash<uint32_t> h;
    return h(zaid.zaid());
  }
};

#endif
