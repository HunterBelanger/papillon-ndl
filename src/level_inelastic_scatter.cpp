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
#include <PapillonNDL/level_inelastic_scatter.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <optional>

namespace pndl {

LevelInelasticScatter::LevelInelasticScatter(const ACE& ace, std::size_t i)
    : C1_(), C2_() {
  C1_ = ace.xss(i);
  C2_ = ace.xss(i + 1);
}

LevelInelasticScatter::LevelInelasticScatter(double Q, double AWR)
    : C1_(), C2_() {
  if (AWR <= 0.) {
    std::string mssg = "AWR must be greater than zero.";
    throw PNDLException(mssg);
  }

  C1_ = (AWR + 1.) * std::abs(Q) / AWR;
  double tmp = AWR / (AWR + 1.);
  C2_ = tmp * tmp;
}

double LevelInelasticScatter::sample_energy(
    double E_in, const std::function<double()>&) const {
  return C2_ * (E_in - C1_);
}

std::optional<double> LevelInelasticScatter::pdf(double /*E_in*/,
                                                 double /*E_out*/) const {
  return std::nullopt;
}

}  // namespace pndl
