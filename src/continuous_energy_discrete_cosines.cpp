/*
 * Copyright 2021, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * Ce logiciel est régi par la licence CeCILL soumise au droit français et
 * respectant les principes de diffusion des logiciels libres. Vous pouvez
 * utiliser, modifier et/ou redistribuer ce programme sous les conditions
 * de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA
 * sur le site "http://www.cecill.info".
 *
 * En contrepartie de l'accessibilité au code source et des droits de copie,
 * de modification et de redistribution accordés par cette licence, il n'est
 * offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 * seule une responsabilité restreinte pèse sur l'auteur du programme,  le
 * titulaire des droits patrimoniaux et les concédants successifs.
 *
 * A cet égard  l'attention de l'utilisateur est attirée sur les risques
 * associés au chargement,  à l'utilisation,  à la modification et/ou au
 * développement et à la reproduction du logiciel par l'utilisateur étant
 * donné sa spécificité de logiciel libre, qui peut le rendre complexe à
 * manipuler et qui le réserve donc à des développeurs et des professionnels
 * avertis possédant  des  connaissances  informatiques approfondies.  Les
 * utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
 * logiciel à leurs besoins dans des conditions permettant d'assurer la
 * sécurité de leurs systèmes et ou de leurs données et, plus généralement,
 * à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
 *
 * Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
 * pris connaissance de la licence CeCILL, et que vous en avez accepté les
 * termes.
 *
 * */
#include <PapillonNDL/continuous_energy_discrete_cosines.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/region_1d.hpp>

#include <iostream>
#include <cmath>

namespace pndl {

ContinuousEnergyDiscreteCosines::ContinuousEnergyDiscreteCosines(const ACE& ace)
    : AngleEnergy(std::make_shared<Region1D>(std::vector<double>({0., 200.}),
                                             std::vector<double>({1., 1.}),
                                             Interpolation::LinLin)),
      incoming_energy_(),
      tables_(),
      Nmu(0) {
  // Make sure the distributions is discrete cosines and energies
  int32_t nxs_7 = ace.nxs(6);
  if (nxs_7 != 2) {
    std::string mssg =
        "ContinuousEnergyDiscreteCosines::ContinuousEnergyDiscreteCosines: The "
        "provided ACE file does not contain a distribution of this form for "
        "Incoherent Inelastic scattering.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  // Read incident energy grid
  int32_t S = ace.jxs(0) - 1;
  uint32_t Ne = ace.xss<uint32_t>(S);  // Number of grid points
  incoming_energy_ = ace.xss(S + 1, Ne);

  // Make sure incident energy grid is sorted
  if (!std::is_sorted(incoming_energy_.begin(), incoming_energy_.end())) {
    std::string mssg =
        "ContinuousEnergyDiscreteCosines::ContinuousEnergyDiscreteCosines: The "
        "incident energy grid is not sorted.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  // Get the number of outgoing discrete cosines
  Nmu = static_cast<uint32_t>(ace.nxs(2)) - 1;

  // Get the starting index for the locators and sizes
  uint32_t i = ace.jxs(2) - 1;
  // Locators in xss for each incident energy
  std::vector<uint32_t> locs = ace.xss<uint32_t>(i, Ne);
  // Number of outgoing energies for each incident energy
  std::vector<uint32_t> Noes = ace.xss<uint32_t>(i + Ne, Ne);

  // Go through all incident energies
  for (size_t ie = 0; ie < Ne; ie++) {
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
    for (size_t oe = 0; oe < Noe; oe++) {
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
      for (size_t m = 0; m < Nmu; m++) {
        tables_.back().cosines[oe][m] = ace.xss(l);
        l++;
      }

      // Check all cosines for the current outgoing energy
      /* The cosines SHOULD all be sorted, but this sometimes isn't the case
         in the ACE files, so I have commented out this check for now. Not
         sure what to do about it for now.
      if(!std::is_sorted(tables_.back().cosines[oe].begin(),
      tables_.back().cosines[oe].end())) { std::string mssg =
      "ContinuousEnergyDiscreteCosines::ContinuousEnergyDiscreteCosines: Cosines
      are not sorted for incoming energy index " + std::to_string(ie) + ",
      outgoing energy index " + std::to_string(oe) + "."; throw
      PNDLException(mssg, __FILE__, __LINE__);
      }*/

      if (tables_.back().cosines[oe].front() < -1.) {
        std::string mssg =
            "ContinuousEnergyDiscreteCosines::ContinuousEnergyDiscreteCosines: "
            "Lowest scattering cosine is less than -1 for incoming energy "
            "index " +
            std::to_string(ie) + ", outgoing energy index " +
            std::to_string(oe) + ".";
        throw PNDLException(mssg, __FILE__, __LINE__);
      }

      if (tables_.back().cosines[oe].back() > 1.) {
        std::string mssg =
            "ContinuousEnergyDiscreteCosines::ContinuousEnergyDiscreteCosines: "
            "Largest scattering cosine is greater than 1 for incoming energy "
            "index " +
            std::to_string(ie) + ", outgoing energy index " +
            std::to_string(oe) + ".";
        throw PNDLException(mssg, __FILE__, __LINE__);
      }
    }  // For all outgoing energies

    // Check outgoing energy distribution
    if (!std::is_sorted(tables_.back().energy.begin(),
                        tables_.back().energy.end())) {
      std::string mssg =
          "ContinuousEnergyDiscreteCosines::ContinuousEnergyDiscreteCosines: "
          "Outgoing energies are not sorted for incoming energy index " +
          std::to_string(ie) + ".";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }

    if (tables_.back().energy.front() < 0.) {
      std::string mssg =
          "ContinuousEnergyDiscreteCosines::ContinuousEnergyDiscreteCosines: "
          "Negative outgoing energies are present for incoming energy index " +
          std::to_string(ie) + ".";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }

    if (!std::is_sorted(tables_.back().cdf.begin(), tables_.back().cdf.end())) {
      std::string mssg =
          "ContinuousEnergyDiscreteCosines::ContinuousEnergyDiscreteCosines: "
          "The outgoing energy CDF is not sorted for incoming energy index " +
          std::to_string(ie) + ".";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }

    for (size_t z = 0; z < Noe; z++) {
      if (tables_.back().pdf[z] < 0.) {
        std::string mssg =
            "ContinuousEnergyDiscreteCosines::ContinuousEnergyDiscreteCosines: "
            "Negative PDF value found for outgoing energy index " +
            std::to_string(z) + ", for incoming energy index " +
            std::to_string(ie) + ".";
        throw PNDLException(mssg, __FILE__, __LINE__);
      }
    }

  }  // For all incident energies
}

AngleEnergyPacket ContinuousEnergyDiscreteCosines::sample_angle_energy(
    double E_in, std::function<double()> rng) const {
  // First we sample the outgoing energy.
  // Determine the index of the bounding tabulated incoming energies
  size_t l;
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
    l = std::distance(incoming_energy_.begin(), in_E_it) - 1;
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
  size_t i, j;        // Indicies for getting the angle latter.
  double xi = rng();  // Random variable for sampling outgoing energy
  if (rng() > f) {
    i = l;
    E_hat = sample_energy(tables_[l], xi, j);
    E_l_1 = E_i_1;
    E_l_M = E_i_M;
  } else {
    i = l + 1;
    E_hat = sample_energy(tables_[l + 1], xi, j);
    E_l_1 = E_i_1_1;
    E_l_M = E_i_1_M;
  }

  double E_out = Emin + ((E_hat - E_l_1) / (E_l_M - E_l_1)) * (Emax - Emin);

  // Now we can go and sample the scattering cosine. This will be done with
  // smearing. First we sample a random cosine index.
  uint32_t k = Nmu * rng();
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

  if (std::abs(mu) > 1.) mu = std::copysign(1., mu);

  return {mu, E_out};
}

}  // namespace pndl
