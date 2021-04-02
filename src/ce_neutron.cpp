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
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/uncorrelated.hpp>
#include <PapillonNDL/constant.hpp>
#include <PapillonNDL/polynomial_1d.hpp>
#include <PapillonNDL/region_1d.hpp>
#include <PapillonNDL/multi_region_1d.hpp>

namespace pndl {

CENeutron::CENeutron(const ACE& ace)
    : zaid_(ace.zaid()),
      awr_(ace.awr()),
      temperature_(ace.temperature()),
      fissile_(ace.fissile()),
      energy_grid_(ace),
      total_xs_(nullptr),
      disappearance_xs_(nullptr),
      elastic_xs_(nullptr),
      photon_production_xs_(nullptr),
      elastic_angle_(nullptr),
      nu_total_(nullptr),
      nu_prompt_(nullptr),
      nu_delayed_(nullptr),
      delayed_groups_(),
      reactions_() {
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
  }

  // Make elastic AngleDistribution
  elastic_angle_ =
      std::make_shared<AngleDistribution>(ace, ace.xss<int>(ace.LAND()));

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
    read_fission_data(ace);
  } else {
    nu_total_ = std::make_shared<Constant>(0.);
    nu_prompt_ = std::make_shared<Constant>(0.);
    nu_delayed_ = std::make_shared<Constant>(0.);
  }
}

CENeutron::CENeutron(const ACE& ace, const CENeutron& nuclide)
    : zaid_(ace.zaid()),
      awr_(ace.awr()),
      temperature_(ace.temperature()),
      fissile_(ace.fissile()),
      energy_grid_(ace),
      total_xs_(nullptr),
      disappearance_xs_(nullptr),
      elastic_xs_(nullptr),
      photon_production_xs_(nullptr),
      elastic_angle_(nullptr),
      nu_total_(nullptr),
      nu_prompt_(nullptr),
      nu_delayed_(nullptr),
      delayed_groups_(),
      reactions_() {
  // Make sure these are the same nuclide !

  if (zaid_ != nuclide.zaid()) {
    std::string mssg =
        "Nuclide::Nuclide: ZAID of ACE doesn't match ZAID of nuclide.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (awr_ != nuclide.awr()) {
    std::string mssg =
        "Nuclide::Nuclide: AWR of ACE doesn't match AWR of nuclide.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

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
  }

  // Copy elastic AngleDistribution
  elastic_angle_ = nuclide.elastic_angle_;

  // Read all reactions
  uint32_t NMT = ace.nxs(3);
  for (uint32_t indx = 0; indx < NMT; indx++) {
    uint32_t MT = ace.xss<uint32_t>(ace.MTR() + indx);

    if (!nuclide.has_reaction(MT)) {
      std::string mssg = "Nuclide::Nuclide: MT = " + std::to_string(MT) +
                         " is present in ACE, but not in nuclide.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }

    Reaction reac(ace, indx, energy_grid_, nuclide.reaction(MT));
    std::pair<uint32_t, Reaction> tmp_pair =
        std::make_pair(MT, Reaction(ace, indx, energy_grid_));
    reactions_.emplace(tmp_pair);
  }

  // Copy fission data from other nuclide
  nu_total_ = nuclide.nu_total_;
  nu_prompt_ = nuclide.nu_prompt_;
  nu_delayed_ = nuclide.nu_delayed_;
  delayed_groups_ = nuclide.delayed_groups_;
}

const EnergyGrid& CENeutron::energy_grid() const { return energy_grid_; }

std::shared_ptr<CrossSection> CENeutron::total_cross_section() const {
  return total_xs_;
}

std::shared_ptr<CrossSection> CENeutron::elastic_cross_section() const {
  return elastic_xs_;
}

std::shared_ptr<CrossSection> CENeutron::disappearance_cross_section() const {
  return disappearance_xs_;
}

std::shared_ptr<CrossSection> CENeutron::photon_production_cross_section()
    const {
  return photon_production_xs_;
}

std::shared_ptr<AngleDistribution> CENeutron::elastic_angle_distribution()
    const {
  return elastic_angle_;
}

void CENeutron::read_fission_data(const ACE& ace) {
  // If prompt and or total neutrons are given
  if(ace.jxs(1) > 0) {
    if (ace.xss(ace.NU()) > 0.) {
      // Either prompt or total given, but not both
      if (ace.DNU() > 0) {  // Prompt is provided, as delayed is present
        nu_prompt_ = read_nu(ace, ace.DNU());
      } else {
        nu_total_ = read_nu(ace, ace.DNU());
      }
    } else {
      // Both prompt and total given
      uint32_t KNU_prmpt = ace.NU() + 1;
      uint32_t KNU_tot = ace.NU() + std::abs(ace.xss<int32_t>(ace.NU())) + 1;

      nu_total_ = read_nu(ace, KNU_tot);
      nu_prompt_ = read_nu(ace, KNU_prmpt);
    }
  }

  // Read delayed nu if given
  if (ace.DNU() > 0) {
    nu_delayed_ = read_nu(ace, ace.DNU());
  }

  // Read all delayed group data
  if (ace.BDD() > 0) {
    uint32_t NGRPS = ace.nxs(7);
    size_t g = 1;
    size_t i = ace.BDD();
    while (g <= NGRPS) {
      delayed_groups_.push_back(DelayedGroup(ace, i, g));
      uint32_t NR = ace.xss<uint32_t>(i + 1);
      uint32_t NE = ace.xss<uint32_t>(i + 2 + 2 * NR);
      i += 3 + 2 * (NR + NE);
      g++;
    }
  }
}

std::shared_ptr<Function1D> CENeutron::read_nu(const ACE& ace, size_t i) {
  uint32_t LNU = ace.xss<uint32_t>(i);

  if (LNU == 1) {  // Polynomial
    return read_polynomial_nu(ace, ++i);
  } else {  // Tabular
    return read_tabular_nu(ace, ++i);
  }
}

std::shared_ptr<Function1D> CENeutron::read_polynomial_nu(const ACE& ace,
                                                            size_t i) {
  uint32_t NC = ace.xss<uint32_t>(i);
  std::vector<double> coeffs = ace.xss(i + 1, NC);
  return std::make_shared<Polynomial1D>(coeffs);
}

std::shared_ptr<Function1D> CENeutron::read_tabular_nu(const ACE& ace,
                                                         size_t i) {
  uint32_t NR = ace.xss<uint32_t>(i);
  uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
  std::vector<double> energy = ace.xss(i + 2 + 2 * NR, NE);
  std::vector<double> y = ace.xss(i + 2 + 2 * NR + NE, NE);

  if (NR == 0 || NR == 1) {
    Interpolation interp = Interpolation::LinLin;
    if (NR == 1) interp = ace.xss<Interpolation>(i + 2);

    return std::make_shared<Region1D>(energy, y, interp);
  } else {
    std::vector<uint32_t> breaks = ace.xss<uint32_t>(i + 1, NR);
    std::vector<Interpolation> interps = ace.xss<Interpolation>(i + 1 + NR, NR);

    return std::make_shared<MultiRegion1D>(breaks, interps, energy, y);
  }
}

}  // namespace pndl
