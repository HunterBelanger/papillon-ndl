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
#ifndef PAPILLON_NDL_DISCRETE_PHOTON_H
#define PAPILLON_NDL_DISCRETE_PHOTON_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/energy_law.hpp>
#include <PapillonNDL/pndl_exception.hpp>

namespace pndl {

/**
 * @brief Energy distribution for discrete photons.
 */
class DiscretePhoton : public EnergyLaw {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in XSS array.
   */
  DiscretePhoton(const ACE& ace, std::size_t i) : lp(), A(), Eg() {
    lp = ace.xss<int>(i);
    Eg = ace.xss(i + 1);
    A = ace.awr();

    if ((lp != 0) && (lp != 1) && (lp != 2)) {
      std::string mssg = "Invalid lp of " + std::to_string(lp) +
                         ". Occurred at index " + std::to_string(i) +
                         " in XSS array.";
      throw PNDLException(mssg);
    }

    if (Eg <= 0.) {
      std::string mssg = "Eg must be greater than zero. Occurred at index " +
                         std::to_string(i) + " in XSS array.";
      throw PNDLException(mssg);
    }

    if (A <= 0.) {
      std::string mssg =
          "Atomic weight ratio must be greater than zero. Occurred at index " +
          std::to_string(i) + " in XSS array.";
      throw PNDLException(mssg);
    }
  }

  /**
   * @param lp Primary indicator flag (0 or 1 is primary, 2 is secondary).
   * @param Eg Energy argument of distribution.
   * @param AWR Atomic Weight Ratio of nuclide.
   */
  DiscretePhoton(int lp, double Eg, double AWR) : lp(lp), A(AWR), Eg(Eg) {
    if ((lp != 0) && (lp != 1) && (lp != 2)) {
      std::string mssg = "Invalid lp of " + std::to_string(lp) + ".";
      throw PNDLException(mssg);
    }

    if (Eg <= 0.) {
      std::string mssg = "Eg must be greater than zero.";
      throw PNDLException(mssg);
    }

    if (A <= 0.) {
      std::string mssg = "Atomic weight ratio must be greater than zero.";
      throw PNDLException(mssg);
    }
  }

  ~DiscretePhoton() = default;

  double sample_energy(double E_in,
                       std::function<double()> /*rng*/) const override final {
    if ((lp == 0) || (lp == 1)) return Eg;
    return Eg + (A / (A + 1.)) * E_in;
  }

  std::optional<double> pdf(double /*E_in*/,
                            double /*E_out*/) const override final {
    return std::nullopt;
  }

  /**
   * @brief Returns the flay indicating whether the photon is a primary or
   * secondary. A secondary photon corresponds with 0 and 1, while a secondary
   * photon has a value of 2.
   */
  int primary_indicator() const { return lp; }

  /**
   * @brief Returns the energy argument for the distribution. If it is a primary
   * photon, this is the outgoing energy, and for a secondary neutron this is
   * the binding energy.
   */
  double photon_energy() const { return Eg; }

 private:
  int lp;
  double A;
  double Eg;
};

}  // namespace pndl

#endif
