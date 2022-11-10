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

#include <ENDFtk.hpp>
#include <boost/hana.hpp>  // Needed for the _c literal for constructing mt4
#include <iostream>
#include <string>
using namespace njoy::ENDFtk;

#include <ndarray.hpp>

#include "constants.hpp"
#include "interpolator.hpp"
#include "tabulated_sab.hpp"

int main(const int argc, const char** argv) {
  std::string fname = "tsl-HinH2O.endf";
  int MAT = 1;

  // std::string fname = "tsl-reactor-graphite-10P.endf";
  // int MAT = 31;

  tree::Tape<std::string> pendf = tree::fromFile(fname);
  file::Type<7> mf7 = pendf.material(MAT).front().file(7).parse<7>();
  section::Type<7, 4> mt4 = mf7.section(4_c);

  auto constants = mt4.constants();

  const int LAT = mt4.LAT();
  const int LASYM = mt4.LASYM();
  const int LLN = constants.LLN();
  const double AWR = constants.AWR()[0];

  std::cout << "File: " << fname << "\n";
  std::cout << "MAT: " << MAT << "\n";
  std::cout << "AWR: " << AWR << "\n";
  std::cout << "LAT: " << LAT << "\n";
  std::cout << "LLN: " << LLN << "\n";
  std::cout << "LASYM: " << LASYM << "\n";

  auto scatteringLaw = mt4.scatteringLaw();

  section::Type<7, 4>::TabulatedFunctions& tsl =
      std::get<section::Type<7, 4>::TabulatedFunctions>(scatteringLaw);

  // Get list of all provided temperatures from the scattering law
  auto lawsForAllTemps = tsl.S();
  std::vector<double> temps = {lawsForAllTemps[0].T().begin(),
                               lawsForAllTemps[0].T().end()};

  section::Type<7, 4>::EffectiveTemperature rawEffectiveTemp =
      mt4.principalEffectiveTemperature();

  // Get Tab1 for the effective temperature
  Tab1 effectiveTemp =
      makeTab1(rawEffectiveTemp.boundaries(), rawEffectiveTemp.interpolants(),
               rawEffectiveTemp.TMOD(), rawEffectiveTemp.TEFF());

  std::vector<TabulatedSab> tsls;
  for (std::size_t i = 0; i < temps.size(); i++) {
    const double T = temps[i];
    const double Teff = effectiveTemp(T);
    tsls.emplace_back(tsl, i, T, Teff, AWR, LAT, LASYM, LLN);
  }

  const double Emin = 1.E-5;
  const double Emax = std::fmax(5., constants.EMAX());

  // Calculate bound xs for a single nuclide of the principal scatterer
  const double xs_b = constants.totalFreeCrossSections()[0] *
                      ((AWR + 1.) / AWR) * ((AWR + 1.) / AWR) /
                      constants.numberAtoms()[0];

  constexpr std::size_t NE = 100;
  std::vector<double> Egrid(NE, 0.);
  const double du =
      (std::log(Emax) - std::log(Emin)) / (static_cast<double>(NE) - 1.);
  for (std::size_t i = 0; i < NE; i++) {
    Egrid[i] = std::exp(std::log(Emin) + static_cast<double>(i) * du);
  }

  if (Egrid.back() > Emax) Egrid.back() = Emax;

  NDArray<double> IIxs({2, NE});
  for (std::size_t i = 0; i < NE; i++) {
    const auto& E = Egrid[i];
    IIxs(0, i) = E;
    const double b_min = Sab::min_beta(E, temps[1]);
    const double b_max = Sab::max_beta(E, temps[1]);
    std::cout << "Running index " << i << ", Energy " << E << " eV.\n";
    IIxs(1, i) = (AWR * xs_b * KB * temps[1] / (4. * E)) *
                 tsls[1].integrate_exp_beta(E, b_min, b_max);
  }

  IIxs.save("iixs.npy");

  return 0;
}
