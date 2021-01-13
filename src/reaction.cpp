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
#include <PapillonNDL/equiprobable_energy_bins.hpp>
#include <PapillonNDL/evaporation.hpp>
#include <PapillonNDL/general_evaporation.hpp>
#include <PapillonNDL/kalbach.hpp>
#include <PapillonNDL/level_inelastic_scatter.hpp>
#include <PapillonNDL/maxwellian.hpp>
#include <PapillonNDL/multi_region_1d.hpp>
#include <PapillonNDL/nbody.hpp>
#include <PapillonNDL/reaction.hpp>
#include <PapillonNDL/region_1d.hpp>
#include <PapillonNDL/tabular_energy.hpp>
#include <PapillonNDL/tabular_energy_angle.hpp>
#include <PapillonNDL/uncorrelated.hpp>
#include <PapillonNDL/watt.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <cmath>
#include <iostream>

namespace pndl {

Reaction::Reaction(const ACE& ace, size_t indx, const EnergyGrid& egrid)
    : mt_(),
      q_(),
      awr_(),
      threshold_(),
      frame_(),
      xs_(),
      angle_energy_(nullptr),
      yield_(nullptr) {
  // Get MT, Q, and AWR
  mt_ = ace.xss<uint32_t>(ace.MTR() + indx);
  q_ = ace.xss(ace.LQR() + indx);
  awr_ = ace.awr();

  // Determine the frame of reference for the outgoing distributions
  frame_ = Frame::Lab;
  if (ace.xss(ace.TYR() + indx) < 0.) frame_ = Frame::CM;

  // Get the yield for the reaction
  double yld = std::abs(ace.xss(ace.TYR() + indx));
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

      yield_ = build_Region1D(energy, y, interp);
    } else {
      std::vector<uint32_t> breaks = ace.xss<uint32_t>(i + 1, NR);
      std::vector<Interpolation> interps =
          ace.xss<Interpolation>(i + 1 + NR, NR);

      yield_ = std::make_shared<MultiRegion1D>(breaks, interps, energy, y);
    }
  }

  // Load the cross section
  uint32_t loca = ace.xss<uint32_t>(ace.LSIG() + indx);
  xs_ = CrossSection(ace, ace.SIG() + loca - 1, egrid);
  threshold_ = xs_.energy(0);

  // Get secondary info if yld != 0 (not an absorption reaction)
  if (yld != 0.) {
    // Get angle distribution location
    int locb = ace.xss<int>(ace.LAND() + indx + 1);

    if (locb >= 0) {
      AngleDistribution angle(ace, locb);

      // Get energy distribution location
      uint32_t locc = ace.xss<uint32_t>(ace.LDLW() + indx);
      size_t i = ace.DLW() + locc - 1;

      // TODO currently ignore extra distribuitons, only read first one
      // Write a warning for now
      if (ace.xss<int>(i) != 0) {
        // there are other distributions
        std::cerr << "\n PapillonNDL WARNING : Reaction MT=" << mt_;
        std::cerr << " for ZAID=" << ace.zaid() << " has multiple";
        std::cerr << " energy distributions.\n";
      }

      int law = ace.xss<int>(i + 1);
      uint32_t idat = ace.xss<uint32_t>(i + 2);
      size_t j = ace.DLW() + idat - 1;

      if (law == 1) {  // Equiprobable Energy Bins
        angle_energy_ = std::make_shared<Uncorrelated>(
            angle, std::make_shared<EquiprobableEnergyBins>(ace, j));

      } else if (law == 3) {  // Level Inelastic Scatter
        angle_energy_ = std::make_shared<Uncorrelated>(
            angle, std::make_shared<LevelInelasticScatter>(ace, j));

      } else if (law == 4) {  // Tabular Energy
        angle_energy_ = std::make_shared<Uncorrelated>(
            angle, std::make_shared<TabularEnergy>(ace, j, ace.DLW()));

      } else if (law == 5) {  // General Evaporation
        angle_energy_ = std::make_shared<Uncorrelated>(
            angle, std::make_shared<GeneralEvaporation>(ace, j));

      } else if (law == 7) {  // Maxwellian
        angle_energy_ = std::make_shared<Uncorrelated>(
            angle, std::make_shared<Maxwellian>(ace, j));

      } else if (law == 9) {  // Evaporation
        angle_energy_ = std::make_shared<Uncorrelated>(
            angle, std::make_shared<Evaporation>(ace, j));

      } else if (law == 11) {  // Watt
        angle_energy_ = std::make_shared<Uncorrelated>(
            angle, std::make_shared<Watt>(ace, j));

      } else if (law == 44) {  // Kalbach
        angle_energy_ = std::make_shared<Kalbach>(ace, j);

      } else if (law == 61) {  // Tabular Energy Angle
        angle_energy_ = std::make_shared<TabularEnergyAngle>(ace, j);

      } else if (law == 66) {  // N-body
        angle_energy_ = std::make_shared<NBody>(ace, j, q_);

      } else {
        // Unknown or unsuported law
        std::string mssg = "Reaction: Unkown energy law " + std::to_string(law);
        mssg += " in reaction MT=" + std::to_string(mt_) + " in ZAID=";
        mssg += std::to_string(ace.zaid());
        throw PNDLException(mssg, __FILE__, __LINE__);
      }
    } else {  // locb < 0, no angle distribution
      // Get energy distribution location
      uint32_t locc = ace.xss<uint32_t>(ace.LDLW() + indx);
      size_t i = ace.DLW() + locc - 1;

      // TODO currently ignore extra distribuitons, only read first one
      // Write a warning for now
      if (ace.xss<int>(i) != 0) {
        // there are other distributions
        std::cerr << "\n PapillonNDL WARNING : Reaction MT=" << mt_;
        std::cerr << " for ZAID=" << ace.zaid() << " has multiple";
        std::cerr << " energy distributions.\n";
      }

      int law = ace.xss<int>(i + 1);
      uint32_t idat = ace.xss<uint32_t>(i + 2);
      size_t j = ace.DLW() + idat - 1;

      if (law == 44) {  // Kalbach
        angle_energy_ = std::make_shared<Kalbach>(ace, j);

      } else if (law == 61) {  // Tabular Energy Angle
        angle_energy_ = std::make_shared<TabularEnergyAngle>(ace, j);

      } else if (law == 66) {  // N-body
        angle_energy_ = std::make_shared<NBody>(ace, j, q_);

      } else {
        // Unknown or unsuported law
        std::string mssg = "Reaction: Unkown energy law " + std::to_string(law);
        mssg += " in reaction MT=" + std::to_string(mt_) + " in ZAID=";
        mssg += std::to_string(ace.zaid());
        throw PNDLException(mssg, __FILE__, __LINE__);
      }
    }
  }
}

const CrossSection& Reaction::cross_section() const { return xs_; }

std::shared_ptr<AngleEnergy> Reaction::angle_energy() const {
  return angle_energy_;
}

std::shared_ptr<Function1D> Reaction::yield() const { return yield_; }

}  // namespace pndl
