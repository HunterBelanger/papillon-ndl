/*
 * Copyright 2020, Hunter Belanger
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
#include <PapillonNDL/tabular_energy_angle.hpp>

#include "constants.hpp"

namespace pndl {

TabularEnergyAngle::TabularEnergyAngle(const ACE& ace, size_t i)
    : incoming_energy_(), tables_() {
  // Get number of interpolation points
  uint32_t NR = ace.xss<uint32_t>(i);
  // Get number of energy points
  uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);

  // Breakpoints and interpolations are not read, as linear-linear
  // interpolation is always used between incoming energies.

  // Read incoming energies
  incoming_energy_ = ace.xss(i + 2 + 2 * NR, NE);
  for (auto& E : incoming_energy_) E *= MEV_TO_EV;

  // Read outgoing energy tables
  for (uint32_t j = 0; j < NE; j++) {
    uint32_t loc = ace.DLW() + ace.xss<uint32_t>(i + 2 + 2 * NR + NE + j) - 1;
    tables_.emplace_back(ace, loc);
  }
}

AngleEnergyPacket TabularEnergyAngle::sample_angle_energy(
    double E_in, std::function<double()> rng) const {
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
    f = interpolation_factor(E_in, incoming_energy_[l],
                             incoming_energy_[l + 1]);
  }

  // Sample outgoing energy, and get R and A for mu
  double E_i_1 = tables_[l].min_energy();
  double E_i_M = tables_[l].max_energy();
  double E_i_1_1 = tables_[l + 1].min_energy();
  double E_i_1_M = tables_[l + 1].max_energy();
  double Emin = E_i_1 + f * (E_i_1_1 - E_i_1);
  double Emax = E_i_M + f * (E_i_1_M - E_i_M);

  AngleEnergyPacket tmp{0., 0.};
  double E_hat = E_in;
  double E_l_1, E_l_M;
  if (rng() > f) {
    tmp = tables_[l].sample_angle_energy(rng);
    E_hat = tmp.energy;
    E_l_1 = E_i_1;
    E_l_M = E_i_M;
  } else {
    tmp = tables_[l + 1].sample_angle_energy(rng);
    E_hat = tmp.energy;
    E_l_1 = E_i_1_1;
    E_l_M = E_i_1_M;
  }

  double E_out = Emin + ((E_hat - E_l_1) / (E_l_M - E_l_1)) * (Emax - Emin);
  double mu = tmp.angle;

  // mu has already ben verified by the EenergyAngleTable to be in interval
  // [-1,1]

  return {mu, E_out};
}

}  // namespace pndl
