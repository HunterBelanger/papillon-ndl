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
#include <PapillonNDL/absorption.hpp>
#include <PapillonNDL/angle_distribution.hpp>
#include <PapillonNDL/cm_distribution.hpp>
#include <PapillonNDL/constant.hpp>
#include <PapillonNDL/discrete_photon.hpp>
#include <PapillonNDL/equiprobable_energy_bins.hpp>
#include <PapillonNDL/evaporation.hpp>
#include <PapillonNDL/general_evaporation.hpp>
#include <PapillonNDL/kalbach.hpp>
#include <PapillonNDL/level_inelastic_scatter.hpp>
#include <PapillonNDL/maxwellian.hpp>
#include <PapillonNDL/multiple_distribution.hpp>
#include <PapillonNDL/nbody.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/reaction_base.hpp>
#include <PapillonNDL/tabular_energy.hpp>
#include <PapillonNDL/tabular_energy_angle.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <PapillonNDL/uncorrelated.hpp>
#include <PapillonNDL/watt.hpp>
#include <cmath>
#include <iostream>
#include <memory>

namespace pndl {

ReactionBase::ReactionBase(const ACE& ace, std::size_t indx)
    : mt_(),
      q_(),
      awr_(),
      threshold_(),
      yield_(nullptr),
      neutron_distribution_(nullptr) {
  // Get MT, Q, and AWR
  mt_ = ace.xss<uint32_t>(ace.MTR() + indx);
  q_ = ace.xss(ace.LQR() + indx);
  awr_ = ace.awr();

  // Determine the frame of reference for the outgoing distributions
  Frame frame_ = Frame::Lab;
  if (ace.xss(ace.TYR() + indx) < 0.) frame_ = Frame::CM;

  // Get the yield for the reaction
  double yld = std::abs(ace.xss(ace.TYR() + indx));
  try {
    if (yld < 100.)
      yield_ = std::make_shared<Constant>(yld);
    else {
      std::size_t i = ace.DLW() + static_cast<uint32_t>(yld) - 101;
      uint32_t NR = ace.xss<uint32_t>(i);
      uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
      std::vector<double> energy = ace.xss(i + 2 + 2 * NR, NE);
      std::vector<double> y = ace.xss(i + 2 + 2 * NR + NE, NE);

      if (NR == 0 || NR == 1) {
        Interpolation interp = Interpolation::LinLin;
        if (NR == 1) interp = ace.xss<Interpolation>(i + 2);

        yield_ = std::make_shared<Tabulated1D>(interp, energy, y);
      } else {
        std::vector<uint32_t> breaks = ace.xss<uint32_t>(i + 1, NR);
        std::vector<Interpolation> interps =
            ace.xss<Interpolation>(i + 1 + NR, NR);

        yield_ = std::make_shared<Tabulated1D>(breaks, interps, energy, y);
      }
    }
  } catch (PNDLException& error) {
    std::string mssg =
        "Could not create yield function for MT = " + std::to_string(mt_) + ".";
    error.add_to_exception(mssg);
    throw error;
  }

  // Here, we set the threshold to zero, just so it has a value, but this should
  // be initalized by the daughter class after initialization.
  threshold_ = 0.;

  // Get secondary info if yld != 0 (not an absorption reaction)
  if (yld != 0.) {
    // Temorary vectors to contain the possbile distributions
    std::vector<std::shared_ptr<AngleEnergy>> distributions;
    std::vector<std::shared_ptr<Tabulated1D>> probabilities;

    try {
      load_neutron_distributions(ace, indx, distributions, probabilities);

      // Make the final distribution
      if (distributions.size() > 1) {
        neutron_distribution_ = std::make_shared<MultipleDistribution>(
            distributions, probabilities);
      } else {
        neutron_distribution_ = distributions.front();
      }

      // Check if we are in the CM frame
      if (frame_ == Frame::CM) {
        std::shared_ptr<AngleEnergy> tmp = neutron_distribution_;
        neutron_distribution_.reset();
        neutron_distribution_ = std::make_shared<CMDistribution>(awr_, q_, tmp);
      }
    } catch (PNDLException& error) {
      std::string mssg =
          "Could not create secondary neutron angle-energy distribution for MT "
          "= " +
          std::to_string(mt_) + ".";
      error.add_to_exception(mssg);
      throw error;
    }
  } else {
    // This is an absorption reaction
    neutron_distribution_ = std::make_shared<Absorption>(mt_);
  }
}

ReactionBase::ReactionBase(uint32_t mt, double q, double awr, double threshold,
                           std::shared_ptr<Function1D> yield,
                           std::shared_ptr<AngleEnergy> neutron_distribution)
    : mt_(mt),
      q_(q),
      awr_(awr),
      threshold_(threshold),
      yield_(yield),
      neutron_distribution_(neutron_distribution) {
  // Make sure the threshold is >= 0
  if (threshold_ < 0.) {
    std::string mssg =
        "Reaction threshold must be greater than or equal to zero.";
    throw PNDLException(mssg);
  }

  if (awr_ <= 0.) {
    std::string mssg = "Atomic weight ratio must be greater than zero.";
    throw PNDLException(mssg);
  }
}

void ReactionBase::load_neutron_distributions(
    const ACE& ace, std::size_t indx,
    std::vector<std::shared_ptr<AngleEnergy>>& distributions,
    std::vector<std::shared_ptr<Tabulated1D>>& probabilities) {
  // Get angle distribution location
  int locb = ace.xss<int>(ace.LAND() + indx + 1);

  std::shared_ptr<AngleDistribution> angle(nullptr);

  // Get angulardistribution if it exists.
  if (locb >= 0) {
    try {
      angle = std::make_shared<AngleDistribution>(ace, locb);
    } catch (PNDLException& error) {
      std::string mssg =
          "Could not create secondary angular distribution for MT = " +
          std::to_string(mt_) + ".";
      error.add_to_exception(mssg);
      throw error;
    }
  }

  // Get energy distribution location
  uint32_t locc = ace.xss<uint32_t>(ace.LDLW() + indx);
  std::size_t i = ace.DLW() + locc - 1;

  // Location of next law (set any non-zero initial value)
  uint32_t lnw = 1;
  while (lnw != 0) {
    lnw = ace.xss<uint32_t>(i);

    // Get law info and location
    int law = ace.xss<int>(i + 1);
    uint32_t idat = ace.xss<uint32_t>(i + 2);
    std::size_t j = ace.DLW() + idat - 1;

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

    probability = std::make_shared<Tabulated1D>(NBT, INT, energy, prob);

    // Build law accordinly
    std::shared_ptr<AngleEnergy> angle_energy(nullptr);
    if (law == 1 && angle) {  // Equiprobable Energy Bins
      angle_energy = std::make_shared<Uncorrelated>(
          *angle, std::make_shared<EquiprobableEnergyBins>(ace, j));

    } else if (law == 2 && angle) {  // Discrete Photon
      angle_energy = std::make_shared<Uncorrelated>(
          *angle, std::make_shared<DiscretePhoton>(ace, j));

    } else if (law == 3 && angle) {  // Level Inelastic Scatter
      angle_energy = std::make_shared<Uncorrelated>(
          *angle, std::make_shared<LevelInelasticScatter>(ace, j));

    } else if (law == 4 && angle) {  // Tabular Energy
      angle_energy = std::make_shared<Uncorrelated>(
          *angle, std::make_shared<TabularEnergy>(ace, j, ace.DLW()));

    } else if (law == 5 && angle) {  // General Evaporation
      angle_energy = std::make_shared<Uncorrelated>(
          *angle, std::make_shared<GeneralEvaporation>(ace, j));

    } else if (law == 7 && angle) {  // Maxwellian
      angle_energy = std::make_shared<Uncorrelated>(
          *angle, std::make_shared<Maxwellian>(ace, j));

    } else if (law == 9 && angle) {  // Evaporation
      angle_energy = std::make_shared<Uncorrelated>(
          *angle, std::make_shared<Evaporation>(ace, j));

    } else if (law == 11 && angle) {  // Watt
      angle_energy = std::make_shared<Uncorrelated>(
          *angle, std::make_shared<Watt>(ace, j));

    } else if (law == 44) {  // Kalbach
      angle_energy = std::make_shared<Kalbach>(ace, j);

    } else if (law == 61) {  // Tabular Energy Angle
      angle_energy = std::make_shared<TabularEnergyAngle>(ace, j, ace.DLW());

    } else if (law == 66) {  // N-body
      angle_energy = std::make_shared<NBody>(ace, j, q_);

    } else if ((law > 0 && law < 6) || law == 7 || law == 9 || law == 11) {
      // Didn't have angle distribution
      std::string mssg = "No anglular distribution provided to acompany law " +
                         std::to_string(law) +
                         " in reaction MT=" + std::to_string(mt_) +
                         " in ZAID=" + std::to_string(ace.zaid().zaid()) + ".";
      throw PNDLException(mssg);
    } else {
      // Unknown or unsuported law
      std::string mssg = "Unkown energy law " + std::to_string(law) +
                         " in reaction MT=" + std::to_string(mt_) +
                         " in ZAID=" + std::to_string(ace.zaid().zaid()) + ".";
      throw PNDLException(mssg);
    }

    // Save law to list of all distributions
    distributions.push_back(angle_energy);
    probabilities.push_back(probability);

    // Get i for next distribution (if there is one)
    i = ace.DLW() + lnw - 1;
  }  // while lnw != 0
}

}  // namespace pndl
