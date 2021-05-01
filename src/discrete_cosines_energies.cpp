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
#include <PapillonNDL/discrete_cosines_energies.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/region_1d.hpp>
#include <algorithm>
#include <cmath>

namespace pndl {

// DiscreteCosinesEnergies can only be used with STIncoherentInelastic,
// and is the only distribution given, so the probability is always 1.
DiscreteCosinesEnergies::DiscreteCosinesEnergies(const ACE& ace)
    : AngleEnergy(std::make_shared<Region1D>(std::vector<double>({0., 200.}),
                                             std::vector<double>({1., 1.}),
                                             Interpolation::LinLin)),
      incoming_energy_(),
      outgoing_energies_(),
      Noe(0),
      Nmu(0),
      skewed_(false) {
  // Make sure the distributions is discrete cosines and energies
  int32_t nxs_7 = ace.nxs(6);

  if (nxs_7 != 0 && nxs_7 != 1) {
    std::string mssg =
        "DiscreteCosinesEnergies::DiscreteCosinesEnergies: The provided ACE "
        "file does not contain a distribution of this form for Incoherent "
        "Inelastic scattering.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (nxs_7 == 1) {
    skewed_ = true;
  }

  // Read incident energy grid
  int32_t S = ace.jxs(0) - 1;
  uint32_t Ne = ace.xss<uint32_t>(S);  // Number of grid points
  incoming_energy_ = ace.xss(S + 1, Ne);

  // Make sure incident energy grid is sorted
  if (!std::is_sorted(incoming_energy_.begin(), incoming_energy_.end())) {
    std::string mssg =
        "DiscreteCosinesEnergies::DiscreteCosinesEnergies: The incident energy "
        "grid is not sorted.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  // Get the number of outgoing discrete energies
  Noe = static_cast<uint32_t>(ace.nxs(3));
  if (skewed_ && Noe < 6) {
    std::string mssg =
        "DiscreteCosinesEnergies::DiscreteCosinesEnergies: A skewed "
        "distribution must have at least 6 outgoing energies. Only " +
        std::to_string(Noe) + " were provided.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  // Get the number of outgoing discrete cosines
  Nmu = static_cast<uint32_t>(ace.nxs(2)) + 1;

  // Get the starting index for the distribution data
  int32_t i = ace.jxs(2) - 1;

  // Go through all incident energies
  for (size_t ie = 0; ie < Ne; ie++) {
    // Add the empy vector of all discrete energies
    outgoing_energies_.push_back({});

    // Get all outgoing energies for the current incident energy
    for (size_t oe = 0; oe < Noe; oe++) {
      double E_out = ace.xss(i);
      // Check E_out
      if (E_out <= 0.) {
        std::string mssg =
            "DiscreteCosinesEnergies::DiscreteCosinesEnergies: Nevative "
            "outgoing energy found at index " +
            std::to_string(i) + ".";
        throw PNDLException(mssg, __FILE__, __LINE__);
      }
      i++;

      std::vector<double> mu = ace.xss(i, Nmu);
      // Check mu grid
      for (size_t j = 0; j < Nmu; j++) {
        if (mu[j] < -1. || mu[j] > 1.) {
          std::string mssg =
              "DiscreteCosinesEnergies::DiscreteCosinesEnergies: Invalid "
              "cosine value found at index " +
              std::to_string(i + j) + ".";
          throw PNDLException(mssg, __FILE__, __LINE__);
        }
      }
      i += Nmu;

      // Add the pair to the last outgoing_energy
      outgoing_energies_.back().push_back({E_out, mu});
    }  // For all outgoing energies
  }    // For all incident energies
}

AngleEnergyPacket DiscreteCosinesEnergies::sample_angle_energy(
    double E_in, std::function<double()> rng) const {
  uint32_t j = 0;
  uint32_t k = static_cast<uint32_t>(Nmu * rng());
  // Sample j for the outgoing energy point
  if (!skewed_) {
    j = static_cast<uint32_t>(Noe * rng());
  } else {
    // Unit probability so that the sum of the probability of all bins
    // is 1 under the skewed conditions.
    double c = 1. / (10. * (static_cast<double>(Noe) - 3.));

    // Sample random value
    double xi = rng();

    if (xi >= 5 * c && xi < 1. - 5 * c) {
      // Select any one of the inner values, all with equal probability.
      // We put the most probable first to make less comparisons.
      j = static_cast<uint32_t>((Noe - 4) * rng() + 2);
    } else if (xi < c) {
      j = 0;
    } else if (xi >= c && xi < 5 * c) {
      j = 1;
    } else if (xi >= 1. - 5 * c && xi < 1. - c) {
      j = Noe - 2;
    } else {
      j = Noe - 1;
    }
  }

  // Find the incoming energy grid point
  auto Eit =
      std::lower_bound(incoming_energy_.begin(), incoming_energy_.end(), E_in);
  if (Eit == incoming_energy_.begin()) {
    // Below lowest incident energy.
    return {outgoing_energies_.front()[j].cosines[k],
            outgoing_energies_.front()[j].energy};
  } else if (Eit == incoming_energy_.end()) {
    // Above largest incident energy.
    return {outgoing_energies_.back()[j].cosines[k],
            outgoing_energies_.back()[j].energy};
  }
  size_t i = std::distance(incoming_energy_.begin(), Eit) - 1;
  double f = (E_in - incoming_energy_[i]) /
             (incoming_energy_[i + 1] - incoming_energy_[i]);

  double E_i_j = outgoing_energies_[i][j].energy;
  double E_i_1_j = outgoing_energies_[i + 1][j].energy;
  double E = E_i_j + f * (E_i_1_j - E_i_j);

  double mu_i_j_k = outgoing_energies_[i][j].cosines[k];
  double mu_i_1_j_k = outgoing_energies_[i + 1][j].cosines[k];
  double mu = mu_i_j_k + f * (mu_i_1_j_k - mu_i_j_k);

  if (std::abs(mu) > 1.) mu = std::copysign(1., mu);

  return {mu, E};
}

}  // namespace pndl
