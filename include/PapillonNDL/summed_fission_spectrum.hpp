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
#ifndef PAPILLON_NDL_SUMMED_FISSION_SPECTRUM_H
#define PAPILLON_NDL_SUMMED_FISSION_SPECTRUM_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/reaction.hpp>
#include <array>
#include <memory>

namespace pndl {

/**
 * @brief In nuclear data evaluations, MT 18 is the reaction which is designated
 *        for fission, and contains the prompt neutron fission spectrum. MT 18
 *        is actually the sum of four other reactions however: MT 19, MT 20, MT
 *        21, and MT 38. These stand for first chance, second chance, third
 *        chance, and fourth chance fission respectively. If these four are
 *        provided insteadd of MT 18, then the prompt fission spectrum is the
 *        average of these four different spectra. This class handles this niche
 * case.
 */
class SummedFissionSpectrum : public AngleEnergy {
 public:
  /**
   * @param mt19 Pointer to STReaction for first chance fission.
   * @param mt20 Pointer to STReaction for first chance fission.
   * @param mt21 Pointer to STReaction for first chance fission.
   * @param mt38 Pointer to STReaction for first chance fission.
   */
  SummedFissionSpectrum(std::shared_ptr<STReaction> mt19,
                        std::shared_ptr<STReaction> mt20,
                        std::shared_ptr<STReaction> mt21,
                        std::shared_ptr<STReaction> mt38);

  AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const override final;

  std::optional<double> angle_pdf(double E_in, double mu) const override final;

  std::optional<double> pdf(double E_in, double mu,
                            double E_out) const override final;

 private:
  std::array<std::shared_ptr<STReaction>, 4> reactions_;

  void compute_probabilities(std::array<double, 4>& probs, double E_in) const;
};

}  // namespace pndl

#endif
