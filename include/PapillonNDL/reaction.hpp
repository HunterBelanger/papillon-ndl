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
#ifndef PAPILLON_NDL_REACTION_H
#define PAPILLON_NDL_REACTION_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/cross_section.hpp>
#include <PapillonNDL/frame.hpp>
#include <PapillonNDL/function_1d.hpp>
#include <memory>

namespace pndl {

/**
 * @brief Holds the cross section and product distributions for a single MT.
 */
class Reaction {
 public:
  /**
   * @param ace ACE file to take reaction from.
   * @param indx Reaction index in the MT array.
   * @param egrid EnergyGrid for the nuclide.
   */
  Reaction(const ACE& ace, size_t indx, const EnergyGrid& egrid);

  /**
   * @param ace ACE file to take cross section from.
   * @param indx Reaction index in the MT array.
   * @param egrid EnergyGrid for the nuclide.
   * @param reac Reaction object to take distributions from.
   */
  Reaction(const ACE& ace, size_t indx, const EnergyGrid& egrid,
           const Reaction& reac);

  ~Reaction() = default;

  /**
   * @brief Returns the MT of the reaction.
   */
  uint32_t mt() const { return mt_; }

  /**
   * @brief Returns the Q-value of the reaction.
   */
  double q() const { return q_; }

  /**
   * @brief Returns the reaction yield for an incident energy.
   * @param E incident energy in MeV.
   */
  double yield(double E) const { return (*yield_)(E); }

  /**
   * @brief Returns the threshold energy for the reaction.
   */
  double threshold() const { return threshold_; }

  /**
   * @brief Returns the frame of referece of the product
   *        distribution data.
   */
  Frame frame() const { return frame_; }

  /**
   * @brief Returns the reaction cross section for a given energy.
   *        Used bisection search.
   * @param E Energy to evaluate the cross section at.
   */
  double xs(double E) const {
    if (E < threshold_) return 0.;

    return xs_->evaluate(E);
  }

  /**
   * @brief Returns the reaction cross section for a given energy.
   * @param E Energy to evaluate the cross section at.
   * @param i Index for the energy grid.
   */
  double xs(double E, size_t i) const {
    if (E < threshold_) return 0.;

    return xs_->evaluate(E, i);
  }

  /**
   * @brief Samples and angle and energy from the reactions product
   *        distribution.
   * @param E_in Incident energy in MeV.
   * @param rng Random number generation function.
   */
  AngleEnergyPacket sample_angle_energy(double E_in,
                                        std::function<double()> rng) const {
    if (distributions_.size() == 0 || E_in < threshold_) return {0., 0.};

    // First select distribution
    double xi = rng();
    double sum = 0.;
    for (size_t d = 0; d < distributions_.size(); d++) {
      sum += distributions_[d]->probability(E_in);
      if (xi < sum) {
        AngleEnergyPacket out =
            distributions_[d]->sample_angle_energy(E_in, rng);
        if (frame_ == Frame::CM) cm_to_lab(E_in, awr_, out);
        return out;
      }
    }

    // Shouldn't get here, but if we do, use the last distribution
    AngleEnergyPacket out =
        distributions_.back()->sample_angle_energy(E_in, rng);
    if (frame_ == Frame::CM) cm_to_lab(E_in, awr_, out);

    return out;
  }

  /**
   * @brief Returns a pointer to the CrossSection for the reaction.
   */
  std::shared_ptr<CrossSection> cross_section() const;

  /**
   * @brief Returns the vector of pointers to all distributions for the
   *        reaction.
   */
  const std::vector<std::shared_ptr<AngleEnergy>>& distributions() const {
    return distributions_;
  }

  /**
   * @brief Returns the number of distributions for the reaction.
   */
  size_t n_distributions() const { return distributions_.size(); }

  /**
   * @brief Returns the ith distributions for the reaction.
   * @param i Index of distribution to fetch.
   */
  std::shared_ptr<AngleEnergy> distribution(size_t i) const {
    return distributions_[i];
  }

  /**
   * @brief Returns a pointer to the function for the reaction yield.
   */
  std::shared_ptr<Function1D> yield() const;

 private:
  uint32_t mt_;
  double q_;
  double awr_;
  double threshold_;
  Frame frame_;
  std::shared_ptr<CrossSection> xs_;
  std::shared_ptr<Function1D> yield_;
  std::vector<std::shared_ptr<AngleEnergy>> distributions_;
};

}  // namespace pndl

#endif
