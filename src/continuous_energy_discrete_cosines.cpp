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
#include <PapillonNDL/continuous_energy_discrete_cosines.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <cmath>
#include <ios>
#include <iostream>
#include <sstream>

namespace pndl {

ContinuousEnergyDiscreteCosines::ContinuousEnergyDiscreteCosines(
    const ACE& ace, bool unit_based_interpolation)
    : incoming_energy_(),
      tables_(),
      Nmu(0),
      unit_based_interpolation_(unit_based_interpolation) {
  // Make sure the distributions is discrete cosines and energies
  int32_t nxs_7 = ace.nxs(6);
  if (nxs_7 != 2) {
    std::string mssg =
        "The provided ACE file does not contain a distribution of this form "
        "for Incoherent Inelastic scattering.";
    throw PNDLException(mssg);
  }

  // Read incident energy grid
  std::size_t S = static_cast<std::size_t>(ace.jxs(0) - 1);
  uint32_t Ne = ace.xss<uint32_t>(S);  // Number of grid points
  incoming_energy_ = ace.xss(S + 1, Ne);

  // Make sure incident energy grid is sorted
  if (!std::is_sorted(incoming_energy_.begin(), incoming_energy_.end())) {
    std::string mssg = "The incident energy grid is not sorted.";
    throw PNDLException(mssg);
  }

  // Get the number of outgoing discrete cosines
  Nmu = static_cast<uint32_t>(ace.nxs(2)) - 1;

  // Get the starting index for the locators and sizes
  uint32_t i = static_cast<uint32_t>(ace.jxs(2) - 1);
  // Locators in xss for each incident energy
  std::vector<uint32_t> locs = ace.xss<uint32_t>(i, Ne);
  // Number of outgoing energies for each incident energy
  std::vector<uint32_t> Noes = ace.xss<uint32_t>(i + Ne, Ne);

  // Go through all incident energies
  for (std::size_t ie = 0; ie < Ne; ie++) {
    // Get the number of outgoing energies for this incident energy
    uint32_t Noe = Noes[ie];
    uint32_t l = locs[ie];

    // Add a blank table to the tables list
    tables_.push_back({{}, {}, {}, {}});

    // Allocate right size for tables
    tables_.back().energy.resize(Noe, 0.);
    tables_.back().pdf.resize(Noe, 0.);
    tables_.back().cdf.resize(Noe, 0.);
    tables_.back().cosines.resize(Noe, std::vector<double>(Nmu, 0.));

    // Go through all outgoing energies
    for (std::size_t oe = 0; oe < Noe; oe++) {
      tables_.back().energy[oe] = ace.xss(l);
      l++;
      tables_.back().pdf[oe] = ace.xss(l);
      // Check that the PDF shouldn't be 0, as sometimes, the last PDF value
      // in the outgoing energy grid is a VERY small negative number
      // (i.e. -2.0E-22). We can safely set these to zero without worry.
      if (tables_.back().pdf[oe] < 0. &&
          std::abs(tables_.back().pdf[oe]) < 1.E-20) {
        tables_.back().pdf[oe] = 0.;
      }
      l++;
      tables_.back().cdf[oe] = ace.xss(l);
      l++;

      // Get all angles
      for (std::size_t m = 0; m < Nmu; m++) {
        tables_.back().cosines[oe][m] = ace.xss(l);
        l++;
      }

      // For ACE files made with NJOY or FRENDY, the cdf usually starts at a
      // value larger than zero. To make sure problems don't arrise, we need
      // to add an entry that is zero. From chapter 7 of the ENDF manual,
      // a cdf of zero corresponds to an exit energy of zero. From Eq. 7.6
      // in the ENDF manual, the pdf for an exit energy of zero, is zero.
      if (oe == 0 && tables_.back().cdf[0] > 0.) {
        tables_.back().energy.insert(tables_.back().energy.begin(), 0.);
        tables_.back().cdf.insert(tables_.back().cdf.begin(), 0.);
        tables_.back().pdf.insert(tables_.back().pdf.begin(), 0.);
        // We now need to create a set is isotropicly distributed discrete
        // angles. An isotropic distribution is used, as very close to Eout=0,
        // alpha is almost constant, which means that S(a,b) will have very
        // little variation, and should therefore provide an isotropic
        // distribution for energies just above 0.
        std::vector<double> discrete_angles(Nmu, 0.);
        const double dmu = 2. / static_cast<double>(Nmu);
        for (uint32_t i_mu = 0; i_mu < Nmu; i_mu++) {
          if (i_mu == 0) {
            discrete_angles[0] = -1. + 0.5 * dmu;
          } else {
            discrete_angles[i_mu] = discrete_angles[i_mu - 1] + dmu;
          }
        }
        tables_.back().cosines.insert(tables_.back().cosines.begin(),
                                      discrete_angles);

        // Advance Noe and oe, due to the added grid point.
        oe++;
        Noe++;
      }

      // Check all cosines for the current outgoing energy.
      // The cosines SHOULD all be sorted, but this sometimes isn't the case
      // in the ACE files. If it isn't, since the angles are all equiprobable, I
      // just sort them.
      if (!std::is_sorted(tables_.back().cosines[oe].begin(),
                          tables_.back().cosines[oe].end())) {
        std::sort(tables_.back().cosines[oe].begin(),
                  tables_.back().cosines[oe].end());
      }

      // All of the cosines should of course be within the interval [-1,1],
      // but I have found thermal scattering laws which is is blatantly not
      // the case. One such example is the light water evaluation at 294K,
      // where there is an upper limit of 1.225, for the ENDF/B-VII.1 data
      // for MCNP.
      if (tables_.back().cosines[oe].front() < -1.) {
        std::stringstream mssg;
        mssg << "Lowest scattering cosine is less than -1 for incoming energy "
                "index "
             << std::to_string(ie) << ", outgoing energy index "
             << std::to_string(oe) << ".";
        throw PNDLException(mssg.str());
      }

      if (tables_.back().cosines[oe].back() > 1.) {
        std::stringstream mssg;
        mssg << "Largest scattering cosine is greater than 1 for incoming "
                "energy index "
             << std::to_string(ie) << ", outgoing energy index "
             << std::to_string(oe) << ".";
        throw PNDLException(mssg.str());
      }
    }  // For all outgoing energies

    // Check outgoing energy distribution
    if (!std::is_sorted(tables_.back().energy.begin(),
                        tables_.back().energy.end())) {
      std::string mssg =
          "Outgoing energies are not sorted for incoming energy index " +
          std::to_string(ie) + ".";
      throw PNDLException(mssg);
    }

    if (tables_.back().energy.front() < 0.) {
      std::string mssg =
          "Negative outgoing energies are present for incoming energy index " +
          std::to_string(ie) + ".";
      throw PNDLException(mssg);
    }

    if (!std::is_sorted(tables_.back().cdf.begin(), tables_.back().cdf.end())) {
      std::string mssg =
          "The outgoing energy CDF is not sorted for incoming energy index " +
          std::to_string(ie) + ".";
      throw PNDLException(mssg);
    }

    for (std::size_t z = 0; z < Noe; z++) {
      // Negative pdf values in the TSLs seems to be systematic. They are always
      // very small in magnitude (usually -1.E-20). Inorder to actually have TSL
      // data to read and use, I simply set it to zero, so long as it wasn't too
      // large a negative value.
      if (tables_.back().pdf[z] < 0.) {
        if (std::abs(tables_.back().pdf[z]) < 1.E-15) {
          tables_.back().pdf[z] = 0.;
        } else {
          std::stringstream mssg;
          mssg << "Large negative PDF value found for incoming energy index "
               << ie << ", and outgoing energy index " << z << ". PDF value is "
               << std::scientific << tables_.back().pdf[z] << ".";
          throw PNDLException(mssg.str());
        }
      }
    }
  }  // For all incident energies
}

AngleEnergyPacket ContinuousEnergyDiscreteCosines::sample_angle_energy(
    double E_in, const std::function<double()>& rng) const {
  if (!unit_based_interpolation_)
    return sample_without_unit_based_interpolation(E_in, rng);
  else
    return sample_with_unit_based_interpolation(E_in, rng);
}

AngleEnergyPacket
ContinuousEnergyDiscreteCosines::sample_with_unit_based_interpolation(
    double E_in, const std::function<double()>& rng) const {
  // First we sample the outgoing energy.
  // Determine the index of the bounding tabulated incoming energies
  std::size_t l;
  double f;  // Interpolation factor
  auto in_E_it =
      std::lower_bound(incoming_energy_.begin(), incoming_energy_.end(), E_in);
  if (in_E_it == incoming_energy_.begin()) {
    l = 0;
    f = 0.;
  } else if (in_E_it == incoming_energy_.end()) {
    l = incoming_energy_.size() - 2;
    f = 1.;
  } else {
    l = static_cast<std::size_t>(
        std::distance(incoming_energy_.begin(), in_E_it) - 1);
    f = (E_in - incoming_energy_[l]) /
        (incoming_energy_[l + 1] - incoming_energy_[l]);
  }

  // Sample outgoing energy like we do with the tabular law, making sure
  // to scale the energy, to make the kinematic valid.
  double E_i_1 = tables_[l].energy.front();
  double E_i_M = tables_[l].energy.back();
  double E_i_1_1 = tables_[l + 1].energy.front();
  double E_i_1_M = tables_[l + 1].energy.back();
  double Emin = E_i_1 + f * (E_i_1_1 - E_i_1);
  double Emax = E_i_M + f * (E_i_1_M - E_i_M);

  double E_hat = E_in;
  double E_l_1, E_l_M;
  std::size_t i, j;   // Indicies for getting the angle latter.
  double xi = rng();  // Random variable for sampling outgoing energy
  if (rng() > f) {
    i = l;
    E_hat = tables_[l].sample_energy(xi, j);
    E_l_1 = E_i_1;
    E_l_M = E_i_M;
  } else {
    i = l + 1;
    E_hat = tables_[l + 1].sample_energy(xi, j);
    E_l_1 = E_i_1_1;
    E_l_M = E_i_1_M;
  }

  double E_out = Emin + ((E_hat - E_l_1) / (E_l_M - E_l_1)) * (Emax - Emin);

  // Now we can go and sample the scattering cosine. This will be done with
  // smearing. First we sample a random cosine index.
  uint32_t k = static_cast<uint32_t>(Nmu * rng());
  f = (xi - tables_[i].cdf[j]) / (tables_[i].cdf[j + 1] - tables_[i].cdf[j]);

  double mu_prime =
      tables_[i].cosines[j][k] +
      f * (tables_[i].cosines[j + 1][k] - tables_[i].cosines[j][k]);

  double mu_left = -1. - (mu_prime + 1.);
  if (k != 0) {
    mu_left =
        tables_[i].cosines[j][k - 1] +
        f * (tables_[i].cosines[j + 1][k - 1] - tables_[i].cosines[j][k - 1]);
  }

  double mu_right = 1. - (mu_prime - 1.);
  if (k != Nmu - 1) {
    mu_right =
        tables_[i].cosines[j][k + 1] +
        f * (tables_[i].cosines[j + 1][k + 1] - tables_[i].cosines[j][k + 1]);
  }

  // Now we smear
  double mu = mu_prime +
              std::min(mu_prime - mu_left, mu_right - mu_prime) * (rng() - 0.5);

  // This is very important as some ACE data has very bad scattering
  // cosine values !
  if (std::abs(mu) > 1.) mu = std::copysign(1., mu);

  return {mu, E_out};
}

AngleEnergyPacket
ContinuousEnergyDiscreteCosines::sample_without_unit_based_interpolation(
    double E_in, const std::function<double()>& rng) const {
  // First we sample the outgoing energy.
  // Determine the index of the bounding tabulated incoming energies
  std::size_t l;
  double f;  // Interpolation factor
  auto in_E_it =
      std::lower_bound(incoming_energy_.begin(), incoming_energy_.end(), E_in);
  if (in_E_it == incoming_energy_.begin()) {
    l = 0;
    f = 0.;
  } else if (in_E_it == incoming_energy_.end()) {
    l = incoming_energy_.size() - 2;
    f = 1.;
  } else {
    l = static_cast<std::size_t>(
        std::distance(incoming_energy_.begin(), in_E_it) - 1);
    f = (E_in - incoming_energy_[l]) /
        (incoming_energy_[l + 1] - incoming_energy_[l]);
  }

  // Always use closest incident energy data
  if (f > 0.5) l++;

  // Sample outgoing energy
  std::size_t i = l, j = 0;
  double xi = rng();
  double E_out = tables_[l].sample_energy(xi, j);

  // No idea why this is done this way in MCNP, Serpent, or OpenMC. A lot of
  // nuclear data was probably made to work right with this though, so we
  // have it here.
  if (E_out < 0.5 * incoming_energy_[l]) {
    E_out *= 2.0 * E_in / incoming_energy_[l] - 1.0;
  } else {
    E_out += E_in - incoming_energy_[l];
  }

  // Now we can go and sample the scattering cosine. This will be done with
  // smearing. First we sample a random cosine index.
  uint32_t k = static_cast<uint32_t>(Nmu * rng());
  f = (xi - tables_[i].cdf[j]) / (tables_[i].cdf[j + 1] - tables_[i].cdf[j]);

  double mu_prime =
      tables_[i].cosines[j][k] +
      f * (tables_[i].cosines[j + 1][k] - tables_[i].cosines[j][k]);

  double mu_left = -1. - (mu_prime + 1.);
  if (k != 0) {
    mu_left =
        tables_[i].cosines[j][k - 1] +
        f * (tables_[i].cosines[j + 1][k - 1] - tables_[i].cosines[j][k - 1]);
  }

  double mu_right = 1. - (mu_prime - 1.);
  if (k != Nmu - 1) {
    mu_right =
        tables_[i].cosines[j][k + 1] +
        f * (tables_[i].cosines[j + 1][k + 1] - tables_[i].cosines[j][k + 1]);
  }

  // Now we smear
  double mu = mu_prime +
              std::min(mu_prime - mu_left, mu_right - mu_prime) * (rng() - 0.5);

  // This is very important as some ACE data has very bad scattering
  // cosine values !
  if (std::abs(mu) > 1.) mu = std::copysign(1., mu);

  return {mu, E_out};
}

double ContinuousEnergyDiscreteCosines::CEDCTable::sample_energy(
    double xi, std::size_t& j) const {
  double E_out = 0.;
  std::size_t l = 0;
  auto cdf_it = std::lower_bound(cdf.begin(), cdf.end(), xi);
  if (cdf_it == cdf.begin()) {
    l = 0;
  } else if (cdf_it == cdf.end()) {
    l = energy.size() - 2;
  } else {
    l = static_cast<std::size_t>(std::distance(cdf.begin(), cdf_it) - 1);
  }

  // Must account for case where pdf_[l] = pdf_[l+1], which means  that
  // the slope is zero, and m=0. This results in nan for the linear alg.
  // To avoid this, must use histogram for that segment.
  if (pdf[l] == pdf[l + 1]) {
    E_out = energy[l] + ((xi - cdf[l]) / pdf[l]);
  } else {
    double m = (pdf[l + 1] - pdf[l]) / (energy[l + 1] - energy[l]);
    E_out =
        energy[l] +
        (1. / m) *
            (std::sqrt(std::max(0., pdf[l] * pdf[l] + 2. * m * (xi - cdf[l]))) -
             pdf[l]);
  }

  // Set j to be l, so we know which cosines to use latter
  // when sampling the angle.
  j = l;

  return E_out;
}

std::optional<double> ContinuousEnergyDiscreteCosines::angle_pdf(
    double /*E_in*/, double /*mu*/) const {
  return std::nullopt;
}

std::optional<double> ContinuousEnergyDiscreteCosines::pdf(
    double /*E_in*/, double /*mu*/, double /*E_out*/) const {
  return std::nullopt;
}

}  // namespace pndl
