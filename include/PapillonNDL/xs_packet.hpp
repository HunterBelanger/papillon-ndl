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
#ifndef PAPILLON_NDL_XSPACKET_H
#define PAPILLON_NDL_XSPACKET_H

/**
 * @file
 * @author Hunter Belanger
 */

namespace pndl {

/**
 * @brief A struct to hold the set of basic cross sections.
 */
struct XSPacket {
  double total;      /**< Total cross section (MT 1) */
  double elastic;    /**< Elastic cross section (MT 2) */
  double inelastic;  /**< Inelastic cross section (MT 3) */
  double absorption; /**< Absorption cross section (MT 27) */
  double fission;    /**< Fission cross section (MT 18) */
  double capture;    /**< Radiative capture cross section (MT 102) */
  double heating;    /**< Heating number */

  XSPacket& operator+=(const XSPacket& other) {
    this->total += other.total;
    this->elastic += other.elastic;
    this->inelastic += other.inelastic;
    this->absorption += other.absorption;
    this->fission += other.fission;
    this->capture += other.capture;
    this->heating += other.heating;
    return *this;
  }

  XSPacket& operator-=(const XSPacket& other) {
    this->total -= other.total;
    this->elastic -= other.elastic;
    this->inelastic -= other.inelastic;
    this->absorption -= other.absorption;
    this->fission -= other.fission;
    this->capture -= other.capture;
    this->heating -= other.heating;
    return *this;
  }

  XSPacket& operator*=(const double& C) {
    this->total *= C;
    this->elastic *= C;
    this->inelastic *= C;
    this->absorption *= C;
    this->fission *= C;
    this->capture *= C;
    this->heating *= C;
    return *this;
  }

  XSPacket& operator/=(const double& C) {
    double D = 1. / C;
    return this->operator*=(D);
  }

  XSPacket operator+(const XSPacket& other) const {
    XSPacket pkt;
    pkt.total = this->total + other.total;
    pkt.elastic = this->elastic + other.elastic;
    pkt.inelastic = this->inelastic + other.inelastic;
    pkt.absorption = this->absorption + other.absorption;
    pkt.fission = this->fission + other.fission;
    pkt.capture = this->capture + other.capture;
    pkt.heating = this->heating + other.heating;
    return pkt;
  }

  XSPacket operator-(const XSPacket& other) const {
    XSPacket pkt;
    pkt.total = this->total - other.total;
    pkt.elastic = this->elastic - other.elastic;
    pkt.inelastic = this->inelastic - other.inelastic;
    pkt.absorption = this->absorption - other.absorption;
    pkt.fission = this->fission - other.fission;
    pkt.capture = this->capture - other.capture;
    pkt.heating = this->heating - other.heating;
    return pkt;
  }

  XSPacket operator*(const double& C) const {
    XSPacket pkt;
    pkt.total = this->total * C;
    pkt.elastic = this->elastic * C;
    pkt.inelastic = this->inelastic * C;
    pkt.absorption = this->absorption * C;
    pkt.fission = this->fission * C;
    pkt.capture = this->capture * C;
    pkt.heating = this->heating * C;
    return pkt;
  }

  XSPacket operator/(const double& C) const {
    double D = 1. / C;
    return this->operator*(D);
  }

  XSPacket operator+() const { return *this; }

  XSPacket operator-() const {
    XSPacket neg;
    neg.total = -this->total;
    neg.elastic = -this->elastic;
    neg.inelastic = -this->inelastic;
    neg.absorption = -this->absorption;
    neg.fission = -this->absorption;
    neg.capture = -this->capture;
    neg.heating = -this->heating;
    return neg;
  }
};

inline XSPacket operator*(const double& C, const XSPacket& xs) {
  return xs * C;
}

}  // namespace pndl

#endif
