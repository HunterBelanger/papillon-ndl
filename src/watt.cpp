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
#include <PapillonNDL/tabulated_1d.hpp>
#include <PapillonNDL/watt.hpp>
#include <cmath>

#include "constants.hpp"

namespace pndl {

Watt::Watt(const ACE& ace, std::size_t i) : a_(), b_(), restriction_energy_() {
  std::size_t original_i = i;
  uint32_t NR = ace.xss<uint32_t>(i);
  uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
  std::vector<uint32_t> NBT_a;
  std::vector<Interpolation> INT_a;

  if (NR == 0) {
    NBT_a = {NE};
    INT_a = {Interpolation::LinLin};
  } else {
    NBT_a = ace.xss<uint32_t>(i + 1, NR);
    INT_a = ace.xss<Interpolation>(i + 1 + NR, NR);
  }

  // Get energy grid
  std::vector<double> energy_a = ace.xss(i + 2 + 2 * NR, NE);
  std::vector<double> a = ace.xss(i + 2 + 2 * NR + NE, NE);

  // Create Function1D pointer
  try {
    a_ = std::make_shared<Tabulated1D>(NBT_a, INT_a, energy_a, a);
  } catch (PNDLException& error) {
    std::string mssg =
        "Could not construct Tabulated1D for the 'a'. Index in the XSS block "
        "is "
        "i = " +
        std::to_string(original_i) + ".";
    error.add_to_exception(mssg);
    throw error;
  }

  // Reset i for b
  i = i + 2 + 2 * NR + 2 * NE;

  NR = ace.xss<uint32_t>(i);
  NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
  std::vector<uint32_t> NBT_b;
  std::vector<Interpolation> INT_b;

  if (NR == 0) {
    NBT_b = {NE};
    INT_b = {Interpolation::LinLin};
  } else {
    NBT_b = ace.xss<uint32_t>(i + 1, NR);
    INT_b = ace.xss<Interpolation>(i + 1 + NR, NR);
  }

  // Get energy grid
  std::vector<double> energy_b = ace.xss(i + 2 + 2 * NR, NE);
  std::vector<double> b = ace.xss(i + 2 + 2 * NR + NE, NE);

  // Create Function1D pointer
  try {
    b_ = std::make_shared<Tabulated1D>(NBT_b, INT_b, energy_b, b);
  } catch (PNDLException& error) {
    std::string mssg =
        "Could not construct Tabulated1D for the 'b'. Index in the XSS block "
        "is "
        "i = " +
        std::to_string(original_i) + ".";
    error.add_to_exception(mssg);
    throw error;
  }

  // Get restriction energy
  restriction_energy_ = ace.xss(i + 2 + 2 * NR + 2 * NE);
}

Watt::Watt(std::shared_ptr<Tabulated1D> a, std::shared_ptr<Tabulated1D> b,
           double restriction_energy)
    : a_(a), b_(b), restriction_energy_(restriction_energy) {}

double Watt::sample_energy(double E_in,
                           const std::function<double()>& rng) const {
  double a = (*a_)(E_in);
  double b = (*b_)(E_in);
  double w = 0.;
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

    w = -a * (std::log(xi1) + std::log(xi2) * c * c);

    E_out = w + 0.25 * a * a * b + (2. * rng() - 1.) * std::sqrt(a * a * b * w);

    if (E_out >= 0. && E_out <= (E_in - restriction_energy_)) sampled = true;
  }

  return E_out;
}

std::optional<double> Watt::pdf(double E_in, double E_out) const {
  double du = E_in - restriction_energy_;
  if (E_out < 0. || E_out > du) return 0.;

  double a = (*a_)(E_in);
  double b = (*b_)(E_in);
  double I = 0.5 * std::sqrt(PI * a * a * a * b / 4.) * std::exp(a * b / 4.);
  I *= std::erf(std::sqrt(du / a) - std::sqrt(a * b / 4.)) +
       std::erf(std::sqrt(du / a) + std::sqrt(a * b / 4.));
  I -= a * std::exp(-du / a) * std::sinh(std::sqrt(b * du));
  return (std::exp(-E_out / a) / I) * std::sinh(std::sqrt(b * E_out));
}

}  // namespace pndl
