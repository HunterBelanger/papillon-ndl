/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2022, Hunter Belanger
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

/**
 * @file
 * @author Hunter Belanger
 */

#include "incoherent_inelastic.hpp"
#include "interpolator.hpp"
#include "tabulated_sab.hpp"
#include "constants.hpp"

#include <boost/hana.hpp>  // Needed for the _c literal for constructing mt4

#include <cmath>
#include <exception>
#include <variant>

IncoherentInelastic::IncoherentInelastic(const section::Type<7, 4>& mt4)
    : sab_(),
      sab_temps_(),
      awr_(0.),
      bound_xs_(0.),
      Emin_(0.),
      Emax_(0.) {
  auto constants = mt4.constants();

  const int LAT = mt4.LAT();
  const int LASYM = mt4.LASYM();
  const int LLN = constants.LLN();
  awr_ = constants.AWR()[0];

  auto scatteringLaw = mt4.scatteringLaw();

  if (std::holds_alternative<section::Type<7, 4>::TabulatedFunctions>(scatteringLaw) == false) {
    throw std::runtime_error("IncoherentInelastic::IncoherentInelastic: No tabulated scattering law in ENDF file.");
  }

  section::Type<7, 4>::TabulatedFunctions& tsl = std::get<section::Type<7, 4>::TabulatedFunctions>(scatteringLaw);

  // Get list of all provided temperatures from the scattering law
  auto lawsForAllTemps = tsl.S();
  sab_temps_ = {lawsForAllTemps[0].T().begin(), lawsForAllTemps[0].T().end()};

  // Get Tab1 for the effective temperature
  section::Type<7, 4>::EffectiveTemperature rawEffectiveTemp = mt4.principalEffectiveTemperature();
  Tab1 effectiveTemp = makeTab1(rawEffectiveTemp.boundaries(), rawEffectiveTemp.interpolants(), rawEffectiveTemp.TMOD(), rawEffectiveTemp.TEFF());

  // Load all S(a,b) scattering laws
  for (std::size_t i = 0; i < sab_temps_.size(); i++) {
    const double T = sab_temps_[i];
    const double Teff = effectiveTemp(T);
    std::unique_ptr<Sab> ST = std::make_unique<TabulatedSab>(tsl, i, T, Teff, awr_, LAT, LASYM, LLN);
    sab_.push_back(std::move(ST));
  }

  // Set min and max energy
  Emin_ = 1.E-5;
  Emax_ = std::fmax(5., constants.EMAX());

  // Calculate bound xs for a single nuclide of the principal scatterer
  bound_xs_ = constants.totalFreeCrossSections()[0] * ((awr_ + 1.) / awr_) * ((awr_ + 1.) / awr_) / constants.numberAtoms()[0];
}

double IncoherentInelastic::ddxs(std::size_t Ti, double Ein, double Eout, double mu) const {
  if (mu < -1. || mu > 1.) {
    throw std::runtime_error("IncoherentInelastic::ddxs: mu must be in interval [-1, 1].");
  }

  const double T = sab_temps_[Ti];
  const Sab& S = *sab_[Ti];
  const double b = (Eout - Ein) / (KB * T);
  const double a =
      (Eout + Ein - 2. * mu * std::sqrt(Ein * Eout)) / (awr_ * KB * T);
  return (awr_ * bound_xs_ * KB * T / (4. * Ein)) * std::exp(-0.5 * b) *
         S(a, b);
}

double IncoherentInelastic::xs(std::size_t Ti, double Ein) const {
  const double T = sab_temps_[Ti];
  const Sab& S = *sab_[Ti];
  const double b_min = Sab::min_beta(Ein, T);
  const double b_max = Sab::max_beta(Ein, T);
  return (awr_ * bound_xs_ * KB * T / (4. * Ein)) *
         S.integrate_exp_beta(Ein, b_min, b_max);
}
