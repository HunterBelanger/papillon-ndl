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
#include <PapillonNDL/cross_section.hpp>
#include <PapillonNDL/reaction_base.hpp>
#include <memory>

namespace pndl {

template <typename XSType>
class Reaction {};

/**
 * @brief Holds all information for a reaction at a single temperature.
 */
template <>
class Reaction<CrossSection> : public ReactionBase {
 public:
  /**
   * @param ace ACE file to take reaction from.
   * @param indx Reaction index in the MT array.
   */
  Reaction(const ACE& ace, std::size_t indx, std::shared_ptr<EnergyGrid> egrid);

  /**
   * @param ace ACE file to take cross section from.
   * @param indx Reaction index in the MT array.
   * @param egrid Pointer to the EnergyGrid for the nuclide.
   * @param reac Reaction object to take distributions from.
   */
  Reaction(const ACE& ace, std::size_t indx, std::shared_ptr<EnergyGrid> egrid,
           const Reaction& reac);

  /**
   * @brief Returns the CrossSection for the reaction.
   */
  const CrossSection& xs() const { return *xs_; }

 private:
  std::shared_ptr<CrossSection> xs_;
};

/**
 * @brief Alias for Reaction<CrossSection>, which contains all data for a
 *        reaction MT at a single temperature.
 */
using STReaction = Reaction<CrossSection>;

}  // namespace pndl

#endif
