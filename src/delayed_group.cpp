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
#include <PapillonNDL/delayed_group.hpp>
#include <PapillonNDL/equiprobable_energy_bins.hpp>
#include <PapillonNDL/evaporation.hpp>
#include <PapillonNDL/general_evaporation.hpp>
#include <PapillonNDL/maxwellian.hpp>
#include <PapillonNDL/multi_region_1d.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/tabular_energy.hpp>
#include <PapillonNDL/watt.hpp>
#include <iostream>

#include "constants.hpp"

namespace pndl {

DelayedGroup::DelayedGroup(const ACE& ace, size_t i, size_t g)
    : decay_constant_(), probability_(nullptr), energy_(nullptr) {
  // Read Decay Constant for group
  decay_constant_ = ace.xss(i);
  decay_constant_ *= SHAKE_TO_SEC;
  i++;

  // Read probability function
  uint32_t NR = ace.xss<uint32_t>(i);
  uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
  std::vector<double> energy = ace.xss(i + 2 + 2 * NR, NE);
  std::vector<double> y = ace.xss(i + 2 + 2 * NR + NE, NE);

  if (NR == 0 || NR == 1) {
    Interpolation interp = Interpolation::LinLin;
    if (NR == 1) interp = ace.xss<Interpolation>(i + 2);

    probability_ = std::make_shared<Region1D>(energy, y, interp);
  } else {
    std::vector<uint32_t> breaks = ace.xss<uint32_t>(i + 1, NR);
    std::vector<Interpolation> interps = ace.xss<Interpolation>(i + 1 + NR, NR);

    probability_ = std::make_shared<MultiRegion1D>(breaks, interps, energy, y);
  }

  // Get energy distribution location
  uint32_t locc = ace.xss<uint32_t>(ace.DNEDL() + g - 1);
  size_t l = ace.DNED() + locc - 1;

  // TODO currently ignore extra distribuitons, only read first one
  // Write a warning for now
  if (ace.xss<int>(l) != 0) {
    // there are other distributions
    std::cerr << "\n PapillonNDL WARNING : Delayed group " + std::to_string(g);
    std::cerr << " for ZAID " << ace.zaid() << " has multiple";
    std::cerr << " energy distributions.\n";
  }

  int law = ace.xss<int>(l + 1);
  uint32_t idat = ace.xss<uint32_t>(l + 2);
  size_t j = ace.DNED() + idat - 1;

  if (law == 1) {  // Equiprobable Energy Bins
    energy_ = std::make_shared<EquiprobableEnergyBins>(ace, j);

  } else if (law == 4) {  // Tabular Energy
    energy_ = std::make_shared<TabularEnergy>(ace, j, ace.DNED());

  } else if (law == 5) {  // General Evaporation
    energy_ = std::make_shared<GeneralEvaporation>(ace, j);

  } else if (law == 7) {  // Maxwellian
    energy_ = std::make_shared<Maxwellian>(ace, j);

  } else if (law == 9) {  // Evaporation
    energy_ = std::make_shared<Evaporation>(ace, j);

  } else if (law == 11) {  // Watt
    energy_ = std::make_shared<Watt>(ace, j);

  } else {
    // Unknown or unsuported law
    std::string mssg = "DelayedGroup: Group " + std::to_string(g) +
                       " has unkown energy law " + std::to_string(law) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }
}

}  // namespace pndl
