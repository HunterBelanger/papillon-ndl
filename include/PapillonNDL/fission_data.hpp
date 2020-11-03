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
#ifndef PAPILLON_NDL_FISSION_DATA_H
#define PAPILLON_NDL_FISSION_DATA_H

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/delayed_group.hpp>
#include <PapillonNDL/function_1d.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <memory>

namespace pndl {

class FissionData {
 public:
  FissionData();
  FissionData(const ACE& ace, std::shared_ptr<AngleEnergy> prmpt);
  ~FissionData() = default;

  std::shared_ptr<Function1D> nu_total() const { return nu_total_; }
  std::shared_ptr<Function1D> nu_prompt() const { return nu_prompt_; }
  std::shared_ptr<Function1D> nu_delayed() const { return nu_delayed_; }

  double nu_total(double E) const {
    if (nu_total_) return (*nu_total_)(E);

    return (*nu_prompt_)(E) + (*nu_delayed_)(E);
  }

  double nu_prompt(double E) const {
    if (nu_prompt_) return (*nu_prompt_)(E);

    // If no prompt, that means no delayed data either, so all
    // neutrons are treated as prompt.
    return (*nu_total_)(E);
  }

  double nu_delayed(double E) const {
    if (nu_delayed_)
      return (*nu_delayed_)(E);

    else if (nu_total_ && nu_prompt_)
      return (*nu_total_)(E) - (*nu_prompt_)(E);

    return 0.;
  }

  size_t ngroups() const { return delayed_groups_.size(); }

  const DelayedGroup& delayed_group(size_t i) const {
    return delayed_groups_[i];
  }

  std::shared_ptr<AngleEnergy> prompt_angle_energy() const {
    return prompt_spectrum_;
  }

  AngleEnergyPacket sample_prompt_angle_energy(double E_in, std::function<double()> rng) const {
    return prompt_spectrum_->sample_angle_energy(E_in, rng);
  }

 private:
  std::shared_ptr<Function1D> nu_total_;
  std::shared_ptr<Function1D> nu_prompt_;
  std::shared_ptr<Function1D> nu_delayed_;

  std::shared_ptr<AngleEnergy> prompt_spectrum_;

  std::vector<DelayedGroup> delayed_groups_;

  std::shared_ptr<Function1D> read_nu(const ACE& ace, size_t i);
  std::shared_ptr<Function1D> read_polynomial_nu(const ACE& ace, size_t i);
  std::shared_ptr<Function1D> read_tabular_nu(const ACE& ace, size_t i);
};

}  // namespace pndl

#endif
