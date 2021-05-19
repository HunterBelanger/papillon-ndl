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
#include <PapillonNDL/multiple_distribution.hpp>
#include <PapillonNDL/pndl_exception.hpp>

namespace pndl {

MultipleDistribution::MultipleDistribution(
    const std::vector<std::shared_ptr<AngleEnergy>>& distributions,
    const std::vector<std::shared_ptr<Tabulated1D>>& probabilities)
    : distributions_(distributions), probabilities_(probabilities) {
  if (distributions_.size() != probabilities_.size()) {
    std::string mssg =
        "MultipleDistribution::MultipleDistribution: Different number of "
        "distributions and probabilities provided.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (distributions_.size() <= 1) {
    std::string mssg =
        "MultipleDistribution::MultipleDistribution: At least two "
        "distributions must be provided.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  for (std::size_t i = 0; i < distributions_.size(); i++) {
    if (!distributions_[i]) {
      std::string mssg =
          "MultipleDistribution::MultipleDistribution: Distribution at index " +
          std::to_string(i) + " is nullptr.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    } else if (!probabilities_[i]) {
      std::string mssg =
          "MultipleDistribution::MultipleDistribution: Probability at index " +
          std::to_string(i) + " is nullptr.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }
}

AngleEnergyPacket MultipleDistribution::sample_angle_energy(
    double E_in, std::function<double()> rng) const {
  // First select distribution
  double xi = rng();
  double sum = 0.;
  for (std::size_t d = 0; d < distributions_.size(); d++) {
    sum += (*probabilities_[d])(E_in);
    if (xi < sum) {
      return distributions_[d]->sample_angle_energy(E_in, rng);
    }
  }

  // Shouldn't get here, but if we do, use the last distribution
  return distributions_.back()->sample_angle_energy(E_in, rng);
}
}  // namespace pndl