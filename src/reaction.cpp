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
#include <PapillonNDL/angle_distribution.hpp>
#include <PapillonNDL/constant.hpp>
#include <PapillonNDL/discrete_photon.hpp>
#include <PapillonNDL/equiprobable_energy_bins.hpp>
#include <PapillonNDL/evaporation.hpp>
#include <PapillonNDL/general_evaporation.hpp>
#include <PapillonNDL/kalbach.hpp>
#include <PapillonNDL/level_inelastic_scatter.hpp>
#include <PapillonNDL/maxwellian.hpp>
#include <PapillonNDL/multi_region_1d.hpp>
#include <PapillonNDL/nbody.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/reaction.hpp>
#include <PapillonNDL/region_1d.hpp>
#include <PapillonNDL/tabular_energy.hpp>
#include <PapillonNDL/tabular_energy_angle.hpp>
#include <PapillonNDL/uncorrelated.hpp>
#include <PapillonNDL/watt.hpp>
#include <cmath>
#include <iostream>

namespace pndl {

Reaction::Reaction(const ACE& ace, size_t indx, const EnergyGrid& egrid)
    : mt_(),
      q_(),
      awr_(),
      threshold_(),
      frame_(),
      xs_(nullptr),
      yield_(nullptr),
      distributions_() {
  // Get MT, Q, and AWR
  mt_ = ace.xss<uint32_t>(ace.MTR() + indx);
  q_ = ace.xss(ace.LQR() + indx);
  awr_ = ace.awr();

  // Determine the frame of reference for the outgoing distributions
  frame_ = Frame::Lab;
  if (ace.xss(ace.TYR() + indx) < 0.) frame_ = Frame::CM;

  // Get the yield for the reaction
  double yld = std::abs(ace.xss(ace.TYR() + indx));
  try {
    if (yld < 100.)
      yield_ = std::make_shared<Constant>(yld);
    else {
      size_t i = ace.DLW() + static_cast<uint32_t>(yld) - 101;
      uint32_t NR = ace.xss<uint32_t>(i);
      uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
      std::vector<double> energy = ace.xss(i + 2 + 2 * NR, NE);
      std::vector<double> y = ace.xss(i + 2 + 2 * NR + NE, NE);

      if (NR == 0 || NR == 1) {
        Interpolation interp = Interpolation::LinLin;
        if (NR == 1) interp = ace.xss<Interpolation>(i + 2);

        yield_ = std::make_shared<Region1D>(energy, y, interp);
      } else {
        std::vector<uint32_t> breaks = ace.xss<uint32_t>(i + 1, NR);
        std::vector<Interpolation> interps =
            ace.xss<Interpolation>(i + 1 + NR, NR);

        yield_ = std::make_shared<MultiRegion1D>(breaks, interps, energy, y);
      }
    }
  } catch (PNDLException& error) {
    std::string mssg =
        "Reaction::Reaction: Could not create yield function for MT = " +
        std::to_string(mt_) + ".";
    error.add_to_exception(mssg, __FILE__, __LINE__);
    throw error;
  }

  // Load the cross section
  try {
    uint32_t loca = ace.xss<uint32_t>(ace.LSIG() + indx);
    xs_ = std::make_shared<CrossSection>(ace, ace.SIG() + loca - 1, egrid);
    threshold_ = xs_->energy(0);
  } catch (PNDLException& error) {
    std::string mssg =
        "Reaction::Reaction: Could not create cross section for MT = " +
        std::to_string(mt_) + ".";
    error.add_to_exception(mssg, __FILE__, __LINE__);
    throw error;
  }

  // Get secondary info if yld != 0 (not an absorption reaction)
  if (yld != 0.) {
    // Get angle distribution location
    int locb = ace.xss<int>(ace.LAND() + indx + 1);

    std::shared_ptr<AngleDistribution> angle(nullptr);

    // Get angulardistribution if it exists.
    if (locb >= 0) {
      try {
        angle = std::make_shared<AngleDistribution>(ace, locb);
      } catch (PNDLException& error) {
        std::string mssg =
            "Reaction::Reaction: Could not create secondary angular "
            "distribution for MT = " +
            std::to_string(mt_) + ".";
        error.add_to_exception(mssg, __FILE__, __LINE__);
        throw error;
      }
    }

    try {
      // Get energy distribution location
      uint32_t locc = ace.xss<uint32_t>(ace.LDLW() + indx);
      size_t i = ace.DLW() + locc - 1;

      // Location of next law (set any non-zero initial value)
      uint32_t lnw = 1;
      while (lnw != 0) {
        lnw = ace.xss<uint32_t>(i);

        // Get law info and location
        int law = ace.xss<int>(i + 1);
        uint32_t idat = ace.xss<uint32_t>(i + 2);
        size_t j = ace.DLW() + idat - 1;

        // Get probability for law
        std::shared_ptr<Tabulated1D> probability(nullptr);
        uint32_t NR = ace.xss<uint32_t>(i + 3);
        uint32_t NE = ace.xss<uint32_t>(i + 4 + 2 * NR);
        std::vector<uint32_t> NBT;
        std::vector<Interpolation> INT;

        if (NR == 0) {
          NBT = {NE};
          INT = {Interpolation::LinLin};
        } else {
          NBT = ace.xss<uint32_t>(i + 4, NR);
          INT = ace.xss<Interpolation>(i + 4 + NR, NR);
        }

        std::vector<double> energy = ace.xss(i + 5 + 2 * NR, NE);
        std::vector<double> prob = ace.xss(i + 5 + 2 * NR + NE, NE);

        if (NBT.size() == 1) {
          probability = std::make_shared<Region1D>(energy, prob, INT[0]);
        } else {
          probability = std::make_shared<MultiRegion1D>(NBT, INT, energy, prob);
        }

        // Build law accordinly
        std::shared_ptr<AngleEnergy> angle_energy(nullptr);
        if (law == 1) {  // Equiprobable Energy Bins
          angle_energy = std::make_shared<Uncorrelated>(
              angle, std::make_shared<EquiprobableEnergyBins>(ace, j),
              probability);

        } else if (law == 2) {  // Discrete Photon
          angle_energy = std::make_shared<Uncorrelated>(
              angle, std::make_shared<DiscretePhoton>(ace, j), probability);

        } else if (law == 3) {  // Level Inelastic Scatter
          angle_energy = std::make_shared<Uncorrelated>(
              angle, std::make_shared<LevelInelasticScatter>(ace, j),
              probability);

        } else if (law == 4) {  // Tabular Energy
          angle_energy = std::make_shared<Uncorrelated>(
              angle, std::make_shared<TabularEnergy>(ace, j, ace.DLW()),
              probability);

        } else if (law == 5) {  // General Evaporation
          angle_energy = std::make_shared<Uncorrelated>(
              angle, std::make_shared<GeneralEvaporation>(ace, j), probability);

        } else if (law == 7) {  // Maxwellian
          angle_energy = std::make_shared<Uncorrelated>(
              angle, std::make_shared<Maxwellian>(ace, j), probability);

        } else if (law == 9) {  // Evaporation
          angle_energy = std::make_shared<Uncorrelated>(
              angle, std::make_shared<Evaporation>(ace, j), probability);

        } else if (law == 11) {  // Watt
          angle_energy = std::make_shared<Uncorrelated>(
              angle, std::make_shared<Watt>(ace, j), probability);

        } else if (law == 44) {  // Kalbach
          angle_energy = std::make_shared<Kalbach>(ace, j, probability);

        } else if (law == 61) {  // Tabular Energy Angle
          angle_energy =
              std::make_shared<TabularEnergyAngle>(ace, j, probability);

        } else if (law == 66) {  // N-body
          angle_energy = std::make_shared<NBody>(ace, j, q_, probability);

        } else {
          // Unknown or unsuported law
          std::string mssg = "Reaction::Reaction: Unkown energy law " +
                             std::to_string(law) +
                             " in reaction MT=" + std::to_string(mt_) +
                             " in ZAID=" + std::to_string(ace.zaid()) + ".";
          throw PNDLException(mssg, __FILE__, __LINE__);
        }

        // Save law to list of all distributions
        distributions_.push_back(angle_energy);

        // Get i for next distribution (if there is one)
        i = ace.DLW() + lnw - 1;
      }  // while lnw != 0
    } catch (PNDLException& error) {
      std::string mssg =
          "Reaction::Reaction: Could not create secondary angle-energy "
          "distribution for MT = " +
          std::to_string(mt_) + ".";
      error.add_to_exception(mssg, __FILE__, __LINE__);
      throw error;
    }
  }
}

Reaction::Reaction(const ACE& ace, size_t indx, const EnergyGrid& egrid,
                   const Reaction& reac)
    : mt_(),
      q_(),
      awr_(),
      threshold_(),
      frame_(),
      xs_(nullptr),
      yield_(nullptr),
      distributions_() {
  // Get MT, Q, and AWR
  mt_ = ace.xss<uint32_t>(ace.MTR() + indx);
  q_ = ace.xss(ace.LQR() + indx);
  awr_ = ace.awr();

  // make sure the MT values agree
  if (mt_ != reac.mt()) {
    std::string mssg =
        "Reaction::Reaction: MT from ACE file doesn't match MT from reaction.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  frame_ = reac.frame();
  yield_ = reac.yield();
  distributions_ = reac.distributions();

  // Get XS from new ACE
  try {
    uint32_t loca = ace.xss<uint32_t>(ace.LSIG() + indx);
    xs_ = std::make_shared<CrossSection>(ace, ace.SIG() + loca - 1, egrid);
    threshold_ = xs_->energy(0);
  } catch (PNDLException& error) {
    std::string mssg =
        "Reaction::Reaction: Could not create cross section for MT = " +
        std::to_string(mt_) + ".";
    error.add_to_exception(mssg, __FILE__, __LINE__);
    throw error;
  }
}

std::shared_ptr<CrossSection> Reaction::cross_section() const { return xs_; }

std::shared_ptr<Function1D> Reaction::yield() const { return yield_; }

}  // namespace pndl
