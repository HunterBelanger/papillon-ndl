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

/**
 * @file
 * @author Hunter Belanger
 */

#include "ace.hpp"

#include <Log.hpp>

#include "constants.hpp"
#include "linearize.hpp"
using namespace njoy;

#include <disco.hpp>
using namespace disco;

#include <array>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <vector>

LinearizedFunction linearize_thermal_scatter_xs(
    const LinearizedFunction& ii, const std::unique_ptr<IncoherentElastic>& ie,
    const std::unique_ptr<CoherentElastic>& ce, const double T) {
  // Initialize an energy grid
  std::vector<double> egrid;

  // Get the number of Bragg Edges, so we can reserve the correct number of
  // points
  const std::size_t NBE = ce ? ce->bragg_edges().size() : 0;

  // Reserve initial size of egrid
  egrid.reserve(ii.x.size() + 2 * NBE);

  // First, add all the Bragg edges, plus the floating point value imediately
  // preceeding each edge. This allows us to pickup the discontinuity.
  if (ce) {
    for (const auto& E : ce->bragg_edges()) {
      egrid.push_back(std::nextafter(E, 0.));
      egrid.push_back(E);
    }
  }

  // Now we can go add all the points from the incoherent inelastic grid.
  for (const auto& E : ii.x) {
    egrid.push_back(E);
  }

  // Sort this grid
  std::sort(egrid.begin(), egrid.end());

  // Now we make a xs object, to evaluate the total xs
  auto txs = [&ii, &ie, &ce, &T](const double& E) {
    double xs_ = ii(E);

    if (ce) xs_ += ce->xs(T, E);

    if (ie) xs_ += ie->xs(T, E);

    return xs_;
  };

  // Create the xs array
  std::vector<double> xs(egrid.size(), 0.);
  for (std::size_t i = 0; i < xs.size(); i++) {
    xs[i] = txs(egrid[i]);
  }

  // Now we linearize and return the result
  return linearize(egrid, xs, txs);
}

void write_ace_ascii(const std::array<std::pair<int32_t, double>, 16>& izaw,
                     const std::array<int32_t, 16>& nxs,
                     const std::array<int32_t, 32>& jxs,
                     const std::vector<double>& xss,
                     const std::vector<bool>& is_float, const std::string& zaid,
                     const double& awr, const double& T,
                     const std::string& comments, const int& mat,
                     const std::filesystem::path& fname) {
  std::string ace;
  auto it = std::back_inserter(ace);

  // Write the header
  using header_line1 = Record<Character<10>,      // ZAID
                              FixedPoint<12, 6>,  // AWR
                              ColumnPosition<1>,  // Blank Space
                              Scientific<11, 4>,  // Temp
                              ColumnPosition<1>,  // Blank Space
                              Character<10>>;     // Processing Date

  using header_line2 = Record<Character<70>,   // Descriptive String
                              Character<10>>;  // Material identifier (mat MAT)

  const double T_MEV = T * K_TO_MEV;

  std::time_t time_t = std::time({});
  std::tm* tm = std::localtime(&time_t);
  std::string date_str;
  date_str.resize(15, ' ');
  std::strftime(&date_str[0], 15, "%d/%m/%Y", tm);
  date_str.resize(10);

  std::stringstream mat_strm;
  mat_strm << "mat " << mat;

  header_line1::write(it, zaid, awr, T_MEV, date_str);
  header_line2::write(it, comments, mat_strm.str());

  // Write the IZAW
  using izaw_line =
      Record<Integer<7>, FixedPoint<11, 0>, Integer<7>, FixedPoint<11, 0>,
             Integer<7>, FixedPoint<11, 0>, Integer<7>, FixedPoint<11, 0>>;

  izaw_line::write(it, izaw[0].first, izaw[0].second, izaw[1].first,
                   izaw[1].second, izaw[2].first, izaw[2].second, izaw[3].first,
                   izaw[3].second);

  izaw_line::write(it, izaw[4].first, izaw[4].second, izaw[5].first,
                   izaw[5].second, izaw[6].first, izaw[6].second, izaw[7].first,
                   izaw[7].second);

  izaw_line::write(it, izaw[8].first, izaw[8].second, izaw[9].first,
                   izaw[9].second, izaw[10].first, izaw[10].second,
                   izaw[11].first, izaw[11].second);

  izaw_line::write(it, izaw[12].first, izaw[12].second, izaw[13].first,
                   izaw[13].second, izaw[14].first, izaw[14].second,
                   izaw[15].first, izaw[15].second);

  // Write the NXS
  using nxs_line = Record<Integer<9>, Integer<9>, Integer<9>, Integer<9>,
                          Integer<9>, Integer<9>, Integer<9>, Integer<9>>;

  nxs_line::write(it, nxs[0], nxs[1], nxs[2], nxs[3], nxs[4], nxs[5], nxs[6],
                  nxs[7]);

  nxs_line::write(it, nxs[8], nxs[9], nxs[10], nxs[11], nxs[12], nxs[13],
                  nxs[14], nxs[15]);

  // Write the JXS
  using jxs_line = Record<Integer<9>, Integer<9>, Integer<9>, Integer<9>,
                          Integer<9>, Integer<9>, Integer<9>, Integer<9>>;

  jxs_line::write(it, jxs[0], jxs[1], jxs[2], jxs[3], jxs[4], jxs[5], jxs[6],
                  jxs[7]);

  jxs_line::write(it, jxs[8], jxs[9], jxs[10], jxs[11], jxs[12], jxs[13],
                  jxs[14], jxs[15]);

  jxs_line::write(it, jxs[16], jxs[17], jxs[18], jxs[19], jxs[20], jxs[21],
                  jxs[22], jxs[23]);

  jxs_line::write(it, jxs[24], jxs[25], jxs[26], jxs[27], jxs[28], jxs[29],
                  jxs[30], jxs[31]);

  // Write the XSS
  for (std::size_t i = 0; i < xss.size(); i++) {
    if (is_float[i]) {
      Scientific<20, 11>::write(xss[i], it);
    } else {
      Integer<20>::write(static_cast<std::int64_t>(xss[i]), it);
    }

    if (i > 0 && (i + 1) % 4 == 0) {
      ace += '\n';
      it = std::back_inserter(ace);
    }
  }

  // Save to file
  std::ofstream file(fname);
  file << ace << std::flush;
}

void write_to_ace(const LinearizedIncoherentInelastic& ii,
                  const std::unique_ptr<IncoherentElastic>& ie,
                  const std::unique_ptr<CoherentElastic>& ce,
                  const std::string& zaid, const double& awr, const double& T,
                  const std::string& comments, const int& mat,
                  const std::filesystem::path& fname) {
  // Initialize data blocks for the ACE file
  std::array<std::pair<int32_t, double>, 16> izaw;
  izaw.fill({0, 0.});

  std::array<int32_t, 16> nxs;
  nxs.fill(0);

  std::array<int32_t, 32> jxs;
  jxs.fill(0);

  std::vector<double> xss;
  std::vector<bool> is_float;

  // First, we will linearize the total thermal scattering xs.
  LinearizedFunction txs =
      linearize_thermal_scatter_xs({ii.egrid, ii.xs}, ie, ce, T);

  // First, we set nxs[4] = 6. This is because nxs(5) is used to indicate one of
  // the following:
  //   * 0 : No Elastic Data
  //   * 3 : Incoherent Elastic Data ONLY
  //   * 4 : Coherent Elastic Data ONLY
  //   * 5 : Mixed Coherent/Incoherent Data
  // We will use the value 6 to indicate that this is an evaluation processed by
  // Panglos, where we use the format outlined bellow in the source.
  //
  // Note that in this version, nxs[1], nxs[2], and nxs[3] are NOT used. nxs[0]
  // holds the total lengths of the xss array, as is standard, and is set at the
  // end of this function.
  nxs[4] = 6;

  //=============================================================================
  // Incoherent Inelastic
  //-----------------------------------------------------------------------------
  // We now write the incoherent inelastic information, starting with the
  // linearized cross section. The start index for the incoherent inelastic
  // energy grid is placed at jxs[0], like in the standard, and the beginning
  // of the xs array in jxs[1]. For this continuous S(a,b) representation, we
  // set nxs[6] = 3;
  xss.reserve(xss.size() + 3 * ii.egrid.size() + 1);
  is_float.reserve(is_float.size() + 3 * ii.egrid.size() + 1);
  jxs[0] = static_cast<int32_t>(xss.size() + 1);
  nxs[6] = 3;

  xss.push_back(
      static_cast<double>(ii.egrid.size()));  // First write number of points
  is_float.push_back(false);
  for (const auto& E : ii.egrid) {
    xss.push_back(E * EV_TO_MEV);  // Now write energy grid
    is_float.push_back(true);
  }

  jxs[1] = static_cast<int32_t>(xss.size() + 1);
  for (const auto& xs : ii.xs) {
    xss.push_back(xs);
    is_float.push_back(true);
  }

  // Starting index of the beta distribution pointers., which we store in
  // jxs[2].
  const std::size_t BptrsStart = xss.size();
  jxs[2] = static_cast<int32_t>(BptrsStart + 1);

  // For each point in the energy grid, we add a "pointer" to the associated
  // beta distribution. We initialize these pointers to all be zero.
  for (std::size_t i = 0; i < ii.egrid.size(); i++) {
    xss.push_back(0.);
    is_float.push_back(false);
  }

  // For each incident energy, we must now write the corresponding beta
  // distribtion, and then the further associated alpha distributions.
  for (std::size_t iE = 0; iE < ii.egrid.size(); iE++) {
    // First, we set the pointer to the distribution, relative to jxs[1] !
    // We +1 for fortran indexing.
    const std::size_t L = xss.size() - static_cast<std::size_t>(jxs[1]) + 1;
    xss[BptrsStart + iE] = static_cast<double>(L);

    // Now we write the number of points in the beta grid.
    xss.push_back(static_cast<double>(ii.beta_dists[iE].beta.size()));
    is_float.push_back(false);

    // Write the beta grid
    for (std::size_t b = 0; b < ii.beta_dists[iE].beta.size(); b++) {
      xss.push_back(ii.beta_dists[iE].beta[b]);
      is_float.push_back(true);
    }
    // Write the pdf
    for (std::size_t b = 0; b < ii.beta_dists[iE].beta.size(); b++) {
      xss.push_back(ii.beta_dists[iE].pdf[b]);
      is_float.push_back(true);
    }
    // Write the cdf
    for (std::size_t b = 0; b < ii.beta_dists[iE].beta.size(); b++) {
      xss.push_back(ii.beta_dists[iE].cdf[b]);
      is_float.push_back(true);
    }

    // Starting index of the Alpa distribution pointers.
    const std::size_t AptrsStart = xss.size();
    // Write the initially empty pointers
    for (std::size_t b = 0; b < ii.beta_dists[iE].beta.size(); b++) {
      xss.push_back(0.);
      is_float.push_back(false);
    }

    // For each beta, we now have an alpha distribution to write
    for (std::size_t iB = 0; iB < ii.beta_dists[iE].beta.size(); iB++) {
      // Write the pointer for the distribution, relative to jxs[1] !
      // We +1 for fortran indexing.
      const std::size_t K = xss.size() - static_cast<std::size_t>(jxs[1]) + 1;
      xss[AptrsStart + iB] = static_cast<double>(K);

      // Now we write the number of points in the alpha grid.
      xss.push_back(
          static_cast<double>(ii.beta_dists[iE].alpha_dists[iB].alpha.size()));
      is_float.push_back(false);

      // Write the alpha grid
      for (std::size_t a = 0;
           a < ii.beta_dists[iE].alpha_dists[iB].alpha.size(); a++) {
        xss.push_back(ii.beta_dists[iE].alpha_dists[iB].alpha[a]);
        is_float.push_back(true);
      }
      // Write the alpha pdf
      for (std::size_t a = 0;
           a < ii.beta_dists[iE].alpha_dists[iB].alpha.size(); a++) {
        xss.push_back(ii.beta_dists[iE].alpha_dists[iB].pdf[a]);
        is_float.push_back(true);
      }
      // Write the alpha cdf
      for (std::size_t a = 0;
           a < ii.beta_dists[iE].alpha_dists[iB].alpha.size(); a++) {
        xss.push_back(ii.beta_dists[iE].alpha_dists[iB].cdf[a]);
        is_float.push_back(true);
      }
    }
  }

  //=============================================================================
  // Coherent Elastic
  //-----------------------------------------------------------------------------
  // If there is Coherent Elastic scattering, the locator goes in jxs[3], like
  // in the standard format. If jxs[3] = 0, then there is no Coherent Elastic
  // scattering
  if (ce) {
    jxs[3] = static_cast<int32_t>(xss.size() + 1);
    const std::size_t NBE = ce->bragg_edges().size();

    // Reserve the needed space for CE
    xss.reserve(xss.size() + 2 * NBE + 1);
    is_float.reserve(is_float.size() + 2 * NBE + 1);

    // First, write number of Bragg Edges
    xss.push_back(static_cast<double>(NBE));
    is_float.push_back(false);

    // Write Bragg edges
    for (const auto& E : ce->bragg_edges()) {
      xss.push_back(E * EV_TO_MEV);
      is_float.push_back(true);
    }

    // Save the index to the structure factor sums in jxs[4].
    jxs[4] = static_cast<int32_t>(xss.size() + 1);

    // Get the structure factors for the desired temperature
    std::vector<double> S = ce->interpolate_structure_factors(T);

    // Write the structure factors
    for (const auto& s : S) {
      xss.push_back(s * EV_TO_MEV);
      is_float.push_back(true);
    }
  } else {
    jxs[3] = 0;
    jxs[4] = 0;
  }

  //=============================================================================
  // Incoherent Elastic
  //-----------------------------------------------------------------------------
  // Place the starting index for Incoherent Elastic at jxs[6]. This aligns with
  // the mixed elastic mode in the new ACE format, but in our custom ACE files,
  // even if only Incoherent Elastic is present, the locator for this channel
  // will be here.
  if (ie) {
    jxs[6] = static_cast<int32_t>(xss.size() + 1);

    // We first write the bound xs, then the Debye-Waller integral for the
    // temperature of interest.
    xss.reserve(xss.size() + 2);
    is_float.reserve(is_float.size() + 2);

    xss.push_back(ie->bound_xs());
    is_float.push_back(true);

    xss.push_back(ie->W()(T));
    is_float.push_back(true);
  } else {
    jxs[6] = 0;
  }

  //=============================================================================
  // TOTAL THERMAL SCATTERING XS
  //-----------------------------------------------------------------------------
  // We now write the total thermal scattering xs to the xss, and record the
  // number of energy points in nxs[8], and the start index in jxs[9]. Neither
  // of these should be in use in the old or new ACE standard for TSLs.
  nxs[8] = static_cast<int32_t>(txs.x.size());
  jxs[9] = static_cast<int32_t>(xss.size() + 1);
  xss.reserve(xss.size() + 2 * txs.x.size());
  is_float.reserve(is_float.size() + 2 * txs.x.size());

  for (const auto& E : txs.x) {
    xss.push_back(E * EV_TO_MEV);
    is_float.push_back(true);
  }

  for (const auto& xs : txs.y) {
    xss.push_back(xs);
    is_float.push_back(true);
  }

  // The total size of the xss goes into nxs[0]
  nxs[0] = static_cast<int32_t>(xss.size());

  // We can now write the ACE file
  write_ace_ascii(izaw, nxs, jxs, xss, is_float, zaid, awr, T, comments, mat,
                  fname);
}
