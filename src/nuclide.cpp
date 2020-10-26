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
#include <PapillonNDL/nuclide.hpp>

namespace pndl {

Nuclide::Nuclide(const ACE& ace)
    : zaid_(ace.zaid()),
      awr_(ace.awr()),
      temperature_(ace.temperature()),
      fissile_(ace.fissile()),
      energy_grid_(ace),
      total_xs_(),
      absorption_xs_(),
      elastic_xs_(),
      elastic_angle_(ace, ace.xss<int>(ace.LAND())),
      reactions_() {
  // Number of energy points
  uint32_t NE = ace.nxs(2);

  total_xs_ = CrossSection(ace, ace.ESZ() + NE, energy_grid_, false);
  absorption_xs_ = CrossSection(ace, ace.ESZ() + 2 * NE, energy_grid_, false);
  elastic_xs_ = CrossSection(ace, ace.ESZ() + 3 * NE, energy_grid_, false);

  // Read all reactions
  uint32_t NMT = ace.nxs(3);
  for (uint32_t indx = 0; indx < NMT; indx++) {
    uint32_t MT = ace.xss<uint32_t>(ace.MTR() + indx);
    Reaction reac(ace, indx, energy_grid_);
    std::pair<uint32_t, Reaction> tmp_pair =
        std::make_pair(MT, Reaction(ace, indx, energy_grid_));
    reactions_.emplace(tmp_pair);
  }

  // TODO read fission data
}

uint32_t Nuclide::ZAID() const { return zaid_; }

double Nuclide::AWR() const { return awr_; }

double Nuclide::temperature() const { return temperature_; }

bool Nuclide::fissile() const { return fissile_; }

const EnergyGrid& Nuclide::energy_grid() const { return energy_grid_; }

const CrossSection& Nuclide::total_cross_section() const { return total_xs_; }

const CrossSection& Nuclide::elastic_cross_section() const {
  return elastic_xs_;
}

const CrossSection& Nuclide::absorption_cross_section() const {
  return absorption_xs_;
}

const AngleDistribution& Nuclide::elastic_angle_distribution() const {
  return elastic_angle_;
}

size_t Nuclide::energy_grid_index(double E) const {
  return energy_grid_.get_lower_index(E);
}

double Nuclide::total_xs(double E) const { return total_xs_(E); }

double Nuclide::total_xs(double E, size_t i) const { return total_xs_(E, i); }

double Nuclide::elastic_xs(double E) const { return elastic_xs_(E); }

double Nuclide::elstic_xs(double E, size_t i) const {
  return elastic_xs_(E, i);
}

double Nuclide::absorption_xs(double E) const { return absorption_xs_(E); }

double Nuclide::absorption_xs(double E, size_t i) const {
  return absorption_xs_(E, i);
}

double Nuclide::sample_elastic_angle(double E,
                                     std::function<double()> rng) const {
  return elastic_angle_.sample_angle(E, rng);
}

bool Nuclide::has_reaction(uint32_t mt) const {
  if (reactions_.find(mt) == reactions_.end()) return false;
  return true;
}

const Reaction& Nuclide::reaction(uint32_t mt) const {
  return reactions_.find(mt)->second;
}

double Nuclide::reaction_xs(uint32_t mt, double E) const {
  if (!has_reaction(mt)) return 0.;

  return reactions_.find(mt)->second.xs(E);
}

double Nuclide::reaction_xs(uint32_t mt, double E, size_t i) const {
  if (!has_reaction(mt)) return 0.;

  return reactions_.find(mt)->second.xs(E, i);
}

}  // namespace pndl
