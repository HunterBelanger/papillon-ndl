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
#include <PapillonNDL/continuous_energy_discrete_cosines.hpp>
#include <PapillonNDL/discrete_cosines_energies.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/st_incoherent_inelastic.hpp>

namespace pndl {

STIncoherentInelastic::STIncoherentInelastic(const ACE& ace,
                                             bool unit_based_interpolation)
    : xs_(nullptr), angle_energy_(nullptr) {
  // Read the XS
  try {
    int32_t S = ace.jxs(0) - 1;
    uint32_t Ne = ace.xss<uint32_t>(S);  // Number of grid points
    std::vector<double> energy = ace.xss(S + 1, Ne);
    std::vector<double> xs = ace.xss(S + 1 + Ne, Ne);
    xs_ = std::make_shared<Region1D>(energy, xs, Interpolation::LinLin);
  } catch (PNDLException& err) {
    std::string mssg = "Could not construct cross section.";
    err.add_to_exception(mssg);
    throw err;
  }

  // Read the angle-energy distribution
  try {
    int32_t nxs_7 = ace.nxs(6);
    if (nxs_7 == 0 || nxs_7 == 1) {
      angle_energy_ = std::make_shared<DiscreteCosinesEnergies>(ace);
    } else if (nxs_7 == 2) {
      angle_energy_ = std::make_shared<ContinuousEnergyDiscreteCosines>(
          ace, unit_based_interpolation);
    } else {
      std::string mssg =
          "Unknown distribution type. Make sure this is a valid thermal "
          "scattering law ACE file.";
      throw PNDLException(mssg);
    }
  } catch (PNDLException& err) {
    std::string mssg = "Could not construct AngleEnergy distribution.";
    err.add_to_exception(mssg);
    throw err;
  }
}
}  // namespace pndl
