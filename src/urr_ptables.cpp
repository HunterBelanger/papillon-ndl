/*
 * Papillon Nuclear Data Library
 * Copyright 2021, Hunter Belanger
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

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/urr_ptables.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <sstream>

namespace pndl {

  URRPTables::URRPTables(const ACE& ace,
      const std::shared_ptr<CrossSection>& elastic,
      const std::shared_ptr<CrossSection>& capture,
      const std::shared_ptr<CrossSection>& fission,
      const std::shared_ptr<CrossSection>& heating,
      const std::vector<STReaction>& reactions):
    interp_(Interpolation::LinLin),
    factors_(false),
    elastic_(nullptr),
    capture_(nullptr),
    fission_(nullptr),
    heating_(nullptr),
    inelastic_(nullptr),
    absorption_(nullptr),
    energy_(nullptr),
    ptables_(nullptr) {
  // Initialize at least the energy_ and ptables_ vectors
  energy_ = std::make_shared<std::vector<double>>();
  ptables_ = std::make_shared<std::vector<PTable>>();

  // Start by checking if there are any ptables or not in the ACE file
  if (ace.jxs(22) == 0) {
    // No URR PTables are given. 
    return;
  }
  
  // Read numbers and flags
  uint32_t Nenergy = ace.xss<uint32_t>(ace.LUNR());
  uint32_t Nptables = ace.xss<uint32_t>(ace.LUNR()+1);
  interp_ = ace.xss<Interpolation>(ace.LUNR()+2);
  int32_t inelastic_flag = ace.xss<int32_t>(ace.LUNR()+3);
  int32_t absorption_flag = ace.xss<int32_t>(ace.LUNR()+4);
  int32_t factors_flag = ace.xss<int32_t>(ace.LUNR()+5);

  // Read energy grid
  *energy_ = ace.xss(ace.LUNR()+6, Nenergy);

  // Read all PTables for each energy
  std::size_t indx = ace.LUNR()+6+Nenergy;
  PTable empty_ptable;
  for (std::size_t ie = 0; ie < Nenergy; ie++) {
    ptables_->push_back(empty_ptable);
    ptables_->back().cdf.resize(Nptables);
    for (std::size_t ipt = 0; ipt < Nptables; ipt++) {
      ptables_->back().cdf[ipt] = ace.xss(indx + ipt); 
    }
    indx += Nptables;

    // Check the CDF for consistency
    if (std::is_sorted(ptables_->back().cdf.begin(), ptables_->back().cdf.end()) == false) {
      std::stringstream mssg;
      mssg << "CDF for incident energy " << ie << " is not sorted."; 
      throw PNDLException(mssg.str());
    }

    if (ptables_->back().cdf[0] <= 0.) {
      std::stringstream mssg;
      mssg << "Initial CDF value for incident energy " << ie << " is less ";
      mssg << "than or equal to zero."; 
      throw PNDLException(mssg.str());
    }

    if (std::abs(1. - ptables_->back().cdf.back()) < 1.E-15) {
      ptables_->back().cdf.back() = 1.;
    } else {
      std::stringstream mssg;
      mssg << "Last CDF value for incident energy " << ie << " is not 1, "; 
      mssg << "but " << ptables_->back().cdf.back() << ".";
      throw PNDLException(mssg.str());
    }
    
    // Get the cross section bands. 
    ptables_->back().xs_bands.resize(Nptables);
    for (std::size_t ipt = 0; ipt < Nptables; ipt++) {
      ptables_->back().xs_bands[ipt].total = ace.xss(indx + ipt); 
    }
    indx += Nptables;

    for (std::size_t ipt = 0; ipt < Nptables; ipt++) {
      ptables_->back().xs_bands[ipt].elastic = ace.xss(indx + ipt); 
    }
    indx += Nptables;

    for (std::size_t ipt = 0; ipt < Nptables; ipt++) {
      ptables_->back().xs_bands[ipt].fission = ace.xss(indx + ipt); 
    }
    indx += Nptables;

    for (std::size_t ipt = 0; ipt < Nptables; ipt++) {
      ptables_->back().xs_bands[ipt].capture = ace.xss(indx + ipt); 
    }
    indx += Nptables;

    for (std::size_t ipt = 0; ipt < Nptables; ipt++) {
      ptables_->back().xs_bands[ipt].heating = ace.xss(indx + ipt); 
    }
    indx += Nptables;
  }

  // Set factors flag
  if (factors_flag == 1) factors_ = true;

  // Check interpolation
  if (interp_ != Interpolation::LinLin &&
      interp_ != Interpolation::LogLog) {
    std::stringstream mssg;
    mssg << "Unsupported interpolation " << interp_ << "." ; 
    throw PNDLException(mssg.str());
  }

  // Set cross sections
  elastic_ = elastic;
  capture_ = capture;
  fission_ = fission;
  heating_ = heating;

  if (inelastic_flag > 0) {
    // Go find relevant MT in reaction list
    for (const auto& reac : reactions) {
      if (reac.mt() == static_cast<uint32_t>(inelastic_flag)) {
        inelastic_ = std::make_shared<CrossSection>(reac.xs()); 
        break;
      }
    }

    if (inelastic_ == nullptr) {
      std::stringstream mssg;
      mssg << "Could not find inelastic MT = " << inelastic_flag << "."; 
      throw PNDLException(mssg.str());
    }
  } /*else if (inelastic_flag == 0) {
    std::string mssg = "Cannot handle inelastic flag of 0."; 
    throw PNDLException(mssg);
  }*/

  if (absorption_flag > 0) {
    // Go find relevant MT in reaction list
    for (const auto& reac : reactions) {
      if (reac.mt() == static_cast<uint32_t>(absorption_flag)) {
        absorption_ = std::make_shared<CrossSection>(reac.xs()); 
        break;
      }
    }

    if (absorption_ == nullptr) {
      std::stringstream mssg;
      mssg << "Could not find absorption MT = " << absorption_flag << "."; 
      throw PNDLException(mssg.str());
    }
  } /*else if (absorption_flag == 0) {
    std::string mssg = "Cannot handle absorption flag of 0."; 
    throw PNDLException(mssg);
  }*/

  /* I can't figure out what to do when inelastic_flag == 0 or when
   * absorption_flag == 0. The ACE manual makes some vague description
   * about calculating the other absoption or inelastic xs using a balance
   * relationship with the smooth cross sections, but I find no other
   * mention of this in any reference. Currently, I do nothing, and I
   * don't account for these extra absorptions or inelastic xs contributions.
   * From what I can tell, this is in agreement with what OpenMC and Scone
   * do, as they both seem to ignore these flags.
   * */
}

}
