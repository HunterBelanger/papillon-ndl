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
#include <PapillonNDL/constant.hpp>
#include <PapillonNDL/fission_data.hpp>
#include <PapillonNDL/multi_region_1d.hpp>
#include <PapillonNDL/polynomial_1d.hpp>
#include <cmath>

namespace pndl {

FissionData::FissionData()
    : awr_(),
      prompt_spectrum_frame_(Frame::Lab),
      nu_total_(nullptr),
      nu_prompt_(nullptr),
      nu_delayed_(nullptr),
      prompt_spectrum_(nullptr),
      delayed_groups_() {}

FissionData::FissionData(const ACE& ace, std::shared_ptr<AngleEnergy> prmpt, Frame frame)
    : awr_(ace.awr()),
      prompt_spectrum_frame_(frame),
      nu_total_(nullptr),
      nu_prompt_(nullptr),
      nu_delayed_(nullptr),
      prompt_spectrum_(prmpt),
      delayed_groups_() {
  if (ace.fissile()) {
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

  } else {
    nu_total_ = std::make_shared<Constant>(0.);
    nu_prompt_ = std::make_shared<Constant>(0.);
    nu_delayed_ = std::make_shared<Constant>(0.);
  }
}

std::shared_ptr<Function1D> FissionData::read_nu(const ACE& ace, size_t i) {
  uint32_t LNU = ace.xss<uint32_t>(i);

  if (LNU == 1) {  // Polynomial
    return read_polynomial_nu(ace, ++i);
  } else {  // Tabular
    return read_tabular_nu(ace, ++i);
  }
}

std::shared_ptr<Function1D> FissionData::read_polynomial_nu(const ACE& ace,
                                                            size_t i) {
  uint32_t NC = ace.xss<uint32_t>(i);
  std::vector<double> coeffs = ace.xss(i + 1, NC);
  return std::make_shared<Polynomial1D>(coeffs);
}

std::shared_ptr<Function1D> FissionData::read_tabular_nu(const ACE& ace,
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
