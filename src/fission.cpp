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

#include <PapillonNDL/absorption.hpp>
#include <PapillonNDL/constant.hpp>
#include <PapillonNDL/difference_1d.hpp>
#include <PapillonNDL/fission.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/polynomial_1d.hpp>
#include <PapillonNDL/sum_1d.hpp>
#include <PapillonNDL/summed_fission_spectrum.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <memory>

namespace pndl {

Fission::Fission(const ACE& ace, std::shared_ptr<EnergyGrid> energy_grid)
    : zaid_(ace.zaid()),
      nu_total_(nullptr),
      nu_prompt_(nullptr),
      nu_delayed_(nullptr),
      mt18_(nullptr),
      mt19_(nullptr),
      mt20_(nullptr),
      mt21_(nullptr),
      mt38_(nullptr),
      prompt_spectrum_(nullptr),
      delayed_families_(),
      mt_list_() {
  if (ace.fissile() == false) {
    nu_total_ = std::make_shared<Constant>(0.);
    nu_prompt_ = std::make_shared<Constant>(0.);
    nu_delayed_ = std::make_shared<Constant>(0.);
    prompt_spectrum_ = std::make_shared<Absorption>(18);
  } else {
    try {
      // Temporary values
      std::shared_ptr<Function1D> total(nullptr);
      std::shared_ptr<Function1D> prompt(nullptr);
      std::shared_ptr<Function1D> delayed(nullptr);

      // If prompt and or total neutrons are given
      if (ace.jxs(1) > 0) {
        if (ace.xss(ace.NU()) > 0.) {
          // Either prompt or total given, but not both
          if (ace.jxs(23) > 0) {
            // Prompt is provided, as delayed is present
            prompt = read_nu(ace, ace.NU());
          } else {
            // No delayed, so this must be total
            total = read_nu(ace, ace.NU());
          }
        } else {
          // Both prompt and total given
          uint32_t KNU_prmpt = ace.NU() + 1;
          uint32_t KNU_tot =
              ace.NU() + std::abs(ace.xss<int32_t>(ace.NU())) + 1;

          total = read_nu(ace, KNU_tot);
          prompt = read_nu(ace, KNU_prmpt);
        }
      }

      // Read delayed nu if given
      if (ace.DNU() > 0) {
        delayed = read_nu(ace, ace.DNU());
      }

      // First we make nu_total_
      if (total) {
        nu_total_ = total;
      } else if (prompt && delayed) {
        nu_total_ = std::make_shared<Sum1D>(prompt, delayed);
      } else if (delayed) {
        nu_total_ = delayed;
      } else {
        // Should never get here, but just in case
        nu_total_ = std::make_shared<Constant>(0.);
      }

      // Now we make nu_prompt_
      if (prompt) {
        nu_prompt_ = prompt;
      } else {
        nu_prompt_ = nu_total_;
      }

      // And finally, delayed
      if (delayed) {
        nu_delayed_ = delayed;
      } else if (total && prompt) {
        nu_delayed_ = std::make_shared<Difference1D>(total, prompt);
      } else {
        nu_delayed_ = std::make_shared<Constant>(0.);
      }
    } catch (PNDLException& err) {
      std::string mssg =
          "Could not construct nu_total, nu_prompt, and nu_delayed.";
      err.add_to_exception(mssg);
      throw err;
    }

    // Read all delayed family data
    try {
      if (ace.BDD() > 0) {
        uint32_t NGRPS = ace.nxs(7);
        std::size_t g = 1;
        std::size_t i = ace.BDD();
        while (g <= NGRPS) {
          delayed_families_.push_back(DelayedFamily(ace, i, g));
          uint32_t NR = ace.xss<uint32_t>(i + 1);
          uint32_t NE = ace.xss<uint32_t>(i + 2 + 2 * NR);
          i += 3 + 2 * (NR + NE);
          g++;
        }
      }
    } catch (PNDLException& err) {
      std::string mssg = "Could not read delayed families.";
      err.add_to_exception(mssg);
      throw err;
    }

    // Read all fission reactions
    try {
      uint32_t NMT = ace.nxs(3);
      mt_list_.reserve(5);
      for (uint32_t indx = 0; indx < NMT; indx++) {
        uint32_t MT = ace.xss<uint32_t>(ace.MTR() + indx);
        if (MT == 18) {
          mt18_ = std::make_shared<STReaction>(ace, indx, energy_grid);
          mt_list_.push_back(18);
        } else if (MT == 19) {
          mt19_ = std::make_shared<STReaction>(ace, indx, energy_grid);
          mt_list_.push_back(19);
        } else if (MT == 20) {
          mt20_ = std::make_shared<STReaction>(ace, indx, energy_grid);
          mt_list_.push_back(20);
        } else if (MT == 21) {
          mt21_ = std::make_shared<STReaction>(ace, indx, energy_grid);
          mt_list_.push_back(21);
        } else if (MT == 38) {
          mt38_ = std::make_shared<STReaction>(ace, indx, energy_grid);
          mt_list_.push_back(38);
        }
      }
    } catch (PNDLException& err) {
      std::string mssg = "Could not read fission reactions.";
      err.add_to_exception(mssg);
      throw err;
    }

    // Create prompt spectrum
    try {
      if (mt18_) {
        prompt_spectrum_ = mt18_->neutron_distribution().shared_from_this();
      } else if (mt19_ || mt20_ || mt21_ || mt38_) {
        prompt_spectrum_ =
            std::make_shared<SummedFissionSpectrum>(mt19_, mt20_, mt21_, mt38_);
      } else {
        // There is no fission apparently. All the nu should have been set to
        // zero, no delayed families should be present, and we will just set the
        // prompt spectrum to absorption.
        prompt_spectrum_ = std::make_shared<Absorption>(18);
      }
    } catch (PNDLException& err) {
      std::string mssg = "Could not create prompt spectrum.";
      err.add_to_exception(mssg);
      throw err;
    }
  }
}

Fission::Fission(const ACE& ace, std::shared_ptr<EnergyGrid> energy_grid,
                 const Fission& fission)
    : zaid_(ace.zaid()),
      nu_total_(fission.nu_total_),
      nu_prompt_(fission.nu_prompt_),
      nu_delayed_(fission.nu_delayed_),
      mt18_(nullptr),
      mt19_(nullptr),
      mt20_(nullptr),
      mt21_(nullptr),
      mt38_(nullptr),
      prompt_spectrum_(nullptr),
      delayed_families_(fission.delayed_families_),
      mt_list_() {
  if (ace.fissile() == false) {
    prompt_spectrum_ = fission.prompt_spectrum_;
  } else {
    // Read all fission reactions
    try {
      uint32_t NMT = ace.nxs(3);
      mt_list_.reserve(5);
      for (uint32_t indx = 0; indx < NMT; indx++) {
        uint32_t MT = ace.xss<uint32_t>(ace.MTR() + indx);
        if (MT == 18) {
          mt18_ = std::make_shared<STReaction>(ace, indx, energy_grid,
                                               *fission.mt18_);
          mt_list_.push_back(18);
        } else if (MT == 19) {
          mt19_ = std::make_shared<STReaction>(ace, indx, energy_grid,
                                               *fission.mt19_);
          mt_list_.push_back(19);
        } else if (MT == 20) {
          mt20_ = std::make_shared<STReaction>(ace, indx, energy_grid,
                                               *fission.mt20_);
          mt_list_.push_back(20);
        } else if (MT == 21) {
          mt21_ = std::make_shared<STReaction>(ace, indx, energy_grid,
                                               *fission.mt21_);
          mt_list_.push_back(21);
        } else if (MT == 38) {
          mt38_ = std::make_shared<STReaction>(ace, indx, energy_grid,
                                               *fission.mt38_);
          mt_list_.push_back(38);
        }
      }
    } catch (PNDLException& err) {
      std::string mssg = "Could not read fission reactions.";
      err.add_to_exception(mssg);
      throw err;
    }

    // Create prompt spectrum
    try {
      if (mt18_) {
        prompt_spectrum_ = mt18_->neutron_distribution().shared_from_this();
      } else if (mt19_ || mt20_ || mt21_ || mt38_) {
        prompt_spectrum_ =
            std::make_shared<SummedFissionSpectrum>(mt19_, mt20_, mt21_, mt38_);
      } else {
        // There is no fission apparently. All the nu should have been set to
        // zero, no delayed families should be present, and we will just set the
        // prompt spectrum to absorption.
        prompt_spectrum_ = std::make_shared<Absorption>(18);
      }
    } catch (PNDLException& err) {
      std::string mssg = "Could not create prompt spectrum.";
      err.add_to_exception(mssg);
      throw err;
    }
  }
}

std::shared_ptr<Function1D> Fission::read_nu(const ACE& ace, std::size_t i) {
  uint32_t LNU = ace.xss<uint32_t>(i);

  if (LNU == 1) {  // Polynomial
    return read_polynomial_nu(ace, ++i);
  } else {  // Tabular
    return read_tabular_nu(ace, ++i);
  }
}

std::shared_ptr<Function1D> Fission::read_polynomial_nu(const ACE& ace,
                                                        std::size_t i) {
  uint32_t NC = ace.xss<uint32_t>(i);
  std::vector<double> coeffs = ace.xss(i + 1, NC);
  return std::make_shared<Polynomial1D>(coeffs);
}

std::shared_ptr<Function1D> Fission::read_tabular_nu(const ACE& ace,
                                                     std::size_t i) {
  uint32_t NR = ace.xss<uint32_t>(i);
  uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
  std::vector<double> energy = ace.xss(i + 2 + 2 * NR, NE);
  std::vector<double> y = ace.xss(i + 2 + 2 * NR + NE, NE);

  if (NR == 0 || NR == 1) {
    Interpolation interp = Interpolation::LinLin;
    if (NR == 1) interp = ace.xss<Interpolation>(i + 2);

    return std::make_shared<Tabulated1D>(interp, energy, y);
  } else {
    std::vector<uint32_t> breaks = ace.xss<uint32_t>(i + 1, NR);
    std::vector<Interpolation> interps = ace.xss<Interpolation>(i + 1 + NR, NR);
    return std::make_shared<Tabulated1D>(breaks, interps, energy, y);
  }
}
}  // namespace pndl
