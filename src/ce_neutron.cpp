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
#include <PapillonNDL/ce_neutron.hpp>

namespace pndl {

CENeutron<CrossSection>::CENeutron(const ACE& ace)
    : CENeutronBase(ace),
      temperature_(ace.temperature()),
      energy_grid_(nullptr),
      total_xs_(nullptr),
      disappearance_xs_(nullptr),
      elastic_xs_(nullptr),
      fission_xs_(nullptr),
      photon_production_xs_(nullptr),
      reactions_() {
  // Construct energy grid
  energy_grid_ = std::make_shared<EnergyGrid>(ace);

  // Number of energy points
  uint32_t NE = ace.nxs(2);
  total_xs_ =
      std::make_shared<CrossSection>(ace, ace.ESZ() + NE, energy_grid_, false);
  disappearance_xs_ = std::make_shared<CrossSection>(ace, ace.ESZ() + 2 * NE,
                                                     energy_grid_, false);
  elastic_xs_ = std::make_shared<CrossSection>(ace, ace.ESZ() + 3 * NE,
                                               energy_grid_, false);

  // Get photon production XS if present
  if (ace.jxs(11) != 0) {
    photon_production_xs_ =
        std::make_shared<CrossSection>(ace, ace.GPD(), energy_grid_, false);
  } else {
    photon_production_xs_ = std::make_shared<CrossSection>(0., energy_grid_);
  }

  // Read all reactions
  uint32_t NMT = ace.nxs(3);
  reaction_indices_.fill(-1);
  int32_t current_reaction_index = 0;
  mt_list_.resize(NMT, 0);
  reactions_.reserve(NMT);
  for (uint32_t indx = 0; indx < mt_list_.size(); indx++) {
    uint32_t MT = ace.xss<uint32_t>(ace.MTR() + indx);
    mt_list_[indx] = MT;
    reactions_.emplace_back(ace, indx, energy_grid_);
    reaction_indices_[MT] = current_reaction_index;
    current_reaction_index++;
  }

  fission_xs_ = compute_fission_xs();
}

CENeutron<CrossSection>::CENeutron(const ACE& ace, const CENeutron& nuclide)
    : CENeutronBase(nuclide),
      temperature_(ace.temperature()),
      energy_grid_(nullptr),
      total_xs_(nullptr),
      disappearance_xs_(nullptr),
      elastic_xs_(nullptr),
      fission_xs_(nullptr),
      photon_production_xs_(nullptr),
      reactions_() {
  // Construct energy grid
  energy_grid_ = std::make_shared<EnergyGrid>(ace);

  // Number of energy points
  uint32_t NE = ace.nxs(2);
  total_xs_ =
      std::make_shared<CrossSection>(ace, ace.ESZ() + NE, energy_grid_, false);
  disappearance_xs_ = std::make_shared<CrossSection>(ace, ace.ESZ() + 2 * NE,
                                                     energy_grid_, false);
  elastic_xs_ = std::make_shared<CrossSection>(ace, ace.ESZ() + 3 * NE,
                                               energy_grid_, false);

  // Get photon production XS if present
  if (ace.jxs(11) != 0) {
    photon_production_xs_ =
        std::make_shared<CrossSection>(ace, ace.GPD(), energy_grid_, false);
  } else {
    photon_production_xs_ = std::make_shared<CrossSection>(0., energy_grid_);
  }

  // Read all reactions
  uint32_t NMT = ace.nxs(3);
  reaction_indices_.fill(-1);
  int32_t current_reaction_index = 0;
  mt_list_.resize(NMT, 0);
  reactions_.reserve(NMT);
  for (uint32_t indx = 0; indx < mt_list_.size(); indx++) {
    uint32_t MT = ace.xss<uint32_t>(ace.MTR() + indx);
    mt_list_[indx] = MT;
    reactions_.emplace_back(ace, indx, energy_grid_);
    reaction_indices_[MT] = current_reaction_index;
    current_reaction_index++;
  }

  fission_xs_ = compute_fission_xs();
}

std::shared_ptr<CrossSection> CENeutron<CrossSection>::compute_fission_xs() {
  if (!fissile_) {
    return std::make_shared<CrossSection>(0., energy_grid_);
  }

  if (this->has_reaction(18)) {
    return std::make_shared<CrossSection>(this->reaction(18).xs());
  }

  // Life is difficult. We need to sum the products.
  std::size_t lowest_index = energy_grid_->size();
  if (this->has_reaction(19)) {
    if (this->reaction(19).xs().index() < lowest_index) {
      lowest_index = this->reaction(19).xs().index();
    }
  }

  if (this->has_reaction(20)) {
    if (this->reaction(20).xs().index() < lowest_index) {
      lowest_index = this->reaction(20).xs().index();
    }
  }

  if (this->has_reaction(21)) {
    if (this->reaction(21).xs().index() < lowest_index) {
      lowest_index = this->reaction(21).xs().index();
    }
  }

  if (this->has_reaction(38)) {
    if (this->reaction(38).xs().index() < lowest_index) {
      lowest_index = this->reaction(38).xs().index();
    }
  }

  // We now know what energy point to start at. We can initialize the vector.
  std::vector<double> fiss_xs(energy_grid_->size() - lowest_index, 0.);

  for (std::size_t i = lowest_index; i < energy_grid_->size(); i++) {
    if (this->has_reaction(19)) {
      fiss_xs[i - lowest_index] += this->reaction(19).xs()[i];
    }

    if (this->has_reaction(20)) {
      fiss_xs[i - lowest_index] += this->reaction(20).xs()[i];
    }

    if (this->has_reaction(21)) {
      fiss_xs[i - lowest_index] += this->reaction(21).xs()[i];
    }

    if (this->has_reaction(38)) {
      fiss_xs[i - lowest_index] += this->reaction(38).xs()[i];
    }
  }

  return std::make_shared<CrossSection>(fiss_xs, energy_grid_, lowest_index);
}

}  // namespace pndl
