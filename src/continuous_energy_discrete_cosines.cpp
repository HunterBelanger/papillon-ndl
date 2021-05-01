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
#include <cmath>
#include <iostream>

namespace pndl {

ContinuousEnergyDiscreteCosines::ContinuousEnergyDiscreteCosines(const ACE& ace, bool unit_based_interpolation)
    : AngleEnergy(std::make_shared<Region1D>(std::vector<double>({0., 200.}),
                                             std::vector<double>({1., 1.}),
                                             Interpolation::LinLin)),
      incoming_energy_(),
      tables_(),
      Nmu(0),
      unit_based_interpolation_(unit_based_interpolation) {
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
      
      /* All of the cosines should of course be within the interval [-1,1],
       * but I have found thermal scattering laws which is is blatantly not
       * the case. One such example is the light water evaluation at 294K,
       * where there is an upper limit of 1.225. Because of this, I have
       * commented out these two checks. To ensure valid results, we must
       * make sure that we validate the cosine range before returning
       * sampled values.
       * */
      /*if (tables_.back().cosines[oe].front() < -1.) {
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
      }*/
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

AngleEnergyPacket ContinuousEnergyDiscreteCosines::sample_angle_energy(double E_in, std::function<double()> rng) const {
  if (!unit_based_interpolation_)
    return sample_without_unit_based_interpolation(E_in, rng);
  else
    return sample_with_unit_based_interpolation(E_in, rng);
}

AngleEnergyPacket ContinuousEnergyDiscreteCosines::sample_with_unit_based_interpolation(
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

AngleEnergyPacket ContinuousEnergyDiscreteCosines::sample_without_unit_based_interpolation(
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
  

  // Always use closest incident energy data
  if (f > 0.5) l++;
  
  // Sample outgoing energy
  size_t i = l, j = 0;
  double xi = rng();
  double E_out = tables_[l].sample_energy(xi, j);
  

  // No idea why this is done this way in MCNP, Serpent, or OpenMC. A lot of 
  // nuclear data was probably made to work right with this though, so we
  // have it here.
  if (E_out < 0.5 * incoming_energy_[l]) {
    E_out *= 2.0*E_in/incoming_energy_[l] - 1.0;
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
    double xi, size_t& j) const {
  double E_out = 0.;
  size_t l = 0;
  auto cdf_it = std::lower_bound(cdf.begin(), cdf.end(), xi);
  if (cdf_it == cdf.begin()) {
    l = 0;
  } else if (cdf_it == cdf.end()) {
    l = energy.size() - 2;
  } else {
    l = std::distance(cdf.begin(), cdf_it) - 1;
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

}  // namespace pndl
