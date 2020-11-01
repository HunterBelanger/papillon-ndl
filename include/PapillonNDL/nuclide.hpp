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
#ifndef PAPILLON_NDL_NUCLIDE_H
#define PAPILLON_NDL_NUCLIDE_H

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_distribution.hpp>
#include <PapillonNDL/fission_data.hpp>
#include <PapillonNDL/reaction.hpp>
#include <unordered_map>

namespace pndl {

class Nuclide {
 public:
  Nuclide(const ACE& ace);
  ~Nuclide() = default;

  uint32_t ZAID() const { return zaid_; }
  double AWR() const { return awr_; }
  double temperature() const { return temperature_; }
  bool fissile() const { return fissile_; }

  const EnergyGrid& energy_grid() const;
  const CrossSection& total_cross_section() const;
  const CrossSection& elastic_cross_section() const;
  const CrossSection& absorption_cross_section() const;
  const AngleDistribution& elastic_angle_distribution() const;

  size_t energy_grid_index(double E) const {
    return energy_grid_.get_lower_index(E);
  }

  double total_xs(double E) const { return total_xs_(E); }

  double total_xs(double E, size_t i) const { return total_xs_(E, i); }

  double elastic_xs(double E) const { return elastic_xs_(E); }

  double elastic_xs(double E, size_t i) const { return elastic_xs_(E, i); }

  double absorption_xs(double E) const { return absorption_xs_(E); }

  double absorption_xs(double E, size_t i) const {
    return absorption_xs_(E, i);
  }

  double sample_elastic_angle(double E, std::function<double()> rng) const {
    return elastic_angle_.sample_angle(E, rng);
  }

  bool has_reaction(uint32_t mt) const {
    if (reactions_.find(mt) == reactions_.end()) return false;
    return true;
  }

  const Reaction& reaction(uint32_t mt) const {
    return reactions_.find(mt)->second;
  }

  double reaction_xs(uint32_t mt, double E) const {
    if (!has_reaction(mt)) return 0.;

    return reactions_.find(mt)->second.xs(E);
  }

  double reaction_xs(uint32_t mt, double E, size_t i) const {
    if (!has_reaction(mt)) return 0.;

    return reactions_.find(mt)->second.xs(E, i);
  }

  const FissionData& fission_data() const { return fission_data_; }

 private:
  uint32_t zaid_;
  double awr_;
  double temperature_;
  bool fissile_;

  EnergyGrid energy_grid_;

  CrossSection total_xs_;
  CrossSection absorption_xs_;
  CrossSection elastic_xs_;

  AngleDistribution elastic_angle_;

  FissionData fission_data_;

  std::unordered_map<uint32_t, Reaction> reactions_;
};

}  // namespace pndl

#endif
