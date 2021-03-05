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
#include <PapillonNDL/nuclide.hpp>
#include <PapillonNDL/uncorrelated.hpp>
#include <PapillonNDL/pndl_exception.hpp>

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
      elastic_angle_(nullptr),
      fission_data_(),
      reactions_() {
  // Number of energy points
  uint32_t NE = ace.nxs(2);

  total_xs_ = CrossSection(ace, ace.ESZ() + NE, energy_grid_, false);
  absorption_xs_ = CrossSection(ace, ace.ESZ() + 2 * NE, energy_grid_, false);
  elastic_xs_ = CrossSection(ace, ace.ESZ() + 3 * NE, energy_grid_, false);

  // Make elastic AngleDistribution
  elastic_angle_ = std::make_shared<AngleDistribution>(ace, ace.xss<int>(ace.LAND()));

  // Read all reactions
  uint32_t NMT = ace.nxs(3);
  for (uint32_t indx = 0; indx < NMT; indx++) {
    uint32_t MT = ace.xss<uint32_t>(ace.MTR() + indx);
    Reaction reac(ace, indx, energy_grid_);
    std::pair<uint32_t, Reaction> tmp_pair =
        std::make_pair(MT, Reaction(ace, indx, energy_grid_));
    reactions_.emplace(tmp_pair);
  }

  if (fissile()) {
    auto Fiss = reaction(18);
    auto prompt_angle_energy = Fiss.angle_energy();
    auto prompt_frame = Fiss.frame();
    fission_data_ = FissionData(ace, prompt_angle_energy, prompt_frame);
  }
}

Nuclide::Nuclide(const ACE& ace, const Nuclide& nuclide)
    : zaid_(ace.zaid()),
      awr_(ace.awr()),
      temperature_(ace.temperature()),
      fissile_(ace.fissile()),
      energy_grid_(ace),
      total_xs_(),
      absorption_xs_(),
      elastic_xs_(),
      elastic_angle_(nullptr),
      fission_data_(),
      reactions_() {

  // Make sure these are the same nuclide !

  if(zaid_ != nuclide.ZAID()) {
    std::string mssg = "Nuclide::Nuclide: ZAID of ACE doesn't match ZAID of nuclide.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if(awr_ != nuclide.AWR()) {
    std::string mssg = "Nuclide::Nuclide: AWR of ACE doesn't match AWR of nuclide.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  // Number of energy points
  uint32_t NE = ace.nxs(2);

  total_xs_ = CrossSection(ace, ace.ESZ() + NE, energy_grid_, false);
  absorption_xs_ = CrossSection(ace, ace.ESZ() + 2 * NE, energy_grid_, false);
  elastic_xs_ = CrossSection(ace, ace.ESZ() + 3 * NE, energy_grid_, false);

  // Copy elastic AngleDistribution
  elastic_angle_ = nuclide.elastic_angle_; 

  // Read all reactions
  uint32_t NMT = ace.nxs(3);
  for (uint32_t indx = 0; indx < NMT; indx++) {
    uint32_t MT = ace.xss<uint32_t>(ace.MTR() + indx);

    if(!nuclide.has_reaction(MT)) {
      std::string mssg = "Nuclide::Nuclide: MT=" + std::to_string(MT);
      mssg += " is present in ACE, but not in nuclide.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }

    Reaction reac(ace, indx, energy_grid_, nuclide.reaction(MT));
    std::pair<uint32_t, Reaction> tmp_pair =
        std::make_pair(MT, Reaction(ace, indx, energy_grid_));
    reactions_.emplace(tmp_pair);
  }

  // Copy fission data from other nuclide
  fission_data_ = nuclide.fission_data_;
}

const EnergyGrid& Nuclide::energy_grid() const { return energy_grid_; }

const CrossSection& Nuclide::total_cross_section() const { return total_xs_; }

const CrossSection& Nuclide::elastic_cross_section() const {
  return elastic_xs_;
}

const CrossSection& Nuclide::absorption_cross_section() const {
  return absorption_xs_;
}

const AngleDistribution& Nuclide::elastic_angle_distribution() const {
  return *elastic_angle_;
}

}  // namespace pndl
