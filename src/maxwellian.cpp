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
#include <PapillonNDL/maxwellian.hpp>
#include <PapillonNDL/multi_region_1d.hpp>
#include <PapillonNDL/region_1d.hpp>
#include <cmath>

#include "constants.hpp"

namespace pndl {

Maxwellian::Maxwellian(const ACE& ace, std::size_t i)
    : temperature_(), restriction_energy_() {
  uint32_t NR = ace.xss<uint32_t>(i);
  uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
  std::vector<uint32_t> NBT;
  std::vector<Interpolation> INT;

  if (NR == 0) {
    NBT = {NE};
    INT = {Interpolation::LinLin};
  } else {
    NBT = ace.xss<uint32_t>(i + 1, NR);
    INT = ace.xss<Interpolation>(i + 1 + NR, NR);
  }

  // Get energy grid
  std::vector<double> energy = ace.xss(i + 2 + 2 * NR, NE);
  std::vector<double> temperature = ace.xss(i + 2 + 2 * NR + NE, NE);

  // Get restriction energy
  restriction_energy_ = ace.xss(i + 2 + 2 * NR + 2 * NE);

  // Create Function1D pointer
  try {
    if (NBT.size() == 1) {
      temperature_ = std::make_unique<Region1D>(energy, temperature, INT[0]);
    } else {
      temperature_ =
          std::make_unique<MultiRegion1D>(NBT, INT, energy, temperature);
    }
  } catch (PNDLException& error) {
    std::string mssg =
        "Could not construct Tabular1D for the effective nuclear temperature. "
        "Index in the XSS block is i = " +
        std::to_string(i) + ".";
    error.add_to_exception(mssg);
    throw error;
  }
}

Maxwellian::Maxwellian(std::shared_ptr<Tabulated1D> temperature,
                       double restriction_energy)
    : temperature_(temperature), restriction_energy_(restriction_energy) {}

double Maxwellian::sample_energy(double E_in,
                                 std::function<double()> rng) const {
  double T = (*temperature_)(E_in);
  double xi1 = 0.;
  double xi2 = 0.;
  double xi3 = 0.;
  double c = 0.;
  double E_out = E_in;
  bool sampled = false;
  while (!sampled) {
    xi1 = rng();
    xi2 = rng();
    xi3 = rng();

    c = std::cos(PI * xi3 / 2.);

    E_out = -T * (std::log(xi1) + std::log(xi2) * c * c);

    if (E_out >= 0. && E_out <= (E_in - restriction_energy_)) sampled = true;
  }

  return E_out;
}

std::optional<double> Maxwellian::pdf(double E_in, double E_out) const {
  double du = E_in - restriction_energy_;
  if (E_out < 0. || E_out > du) return 0.;

  double T = (*temperature_)(E_in);
  double I = std::pow(T, 3. / 2.) * ((std::sqrt(PI) / 2.) * std::erf(du / T) -
                                     std::sqrt(du / T) * std::exp(-du / T));
  return (std::sqrt(E_out) / I) * std::exp(-E_out / T);
}

}  // namespace pndl
