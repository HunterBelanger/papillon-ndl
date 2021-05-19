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
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/st_incoherent_elastic.hpp>

namespace pndl {

STIncoherentElastic::STIncoherentElastic(const ACE& ace)
    : xs_(nullptr), Nmu(0), incoming_energy_(), cosines_() {
  // Fist make sure ACE file does indeed give coherent elastic scattering
  int32_t elastic_mode = ace.nxs(4);
  if (elastic_mode == 4 || ace.jxs(3) == 0) {
    // Make a zero xs incase a user try to get the XS
    std::vector<double> E(2, 0.);
    E[1] = 100.;
    std::vector<double> xs_vals(2, 0.);
    xs_ = std::make_shared<Region1D>(E, xs_vals, Interpolation::Histogram);
  } else {
    // Get index to incident energy
    int32_t i = ace.jxs(3) - 1;
    uint32_t Ne = ace.xss<uint32_t>(i);
    incoming_energy_ = ace.xss(i + 1, Ne);
    std::vector<double> xs_vals = ace.xss(i + 1 + Ne, Ne);

    // Make sure no negative XS values
    for (size_t j = 0; j < xs_vals.size(); j++) {
      if (xs_vals[j] < 0.) {
        std::string mssg =
            "STIncoherentElastic::STIncoherentElastic: Negative cross section "
            "found at index " +
            std::to_string(j) + ".";
        throw PNDLException(mssg, __FILE__, __LINE__);
      }
    }

    xs_ = std::make_shared<Region1D>(incoming_energy_, xs_vals,
                                     Interpolation::LinLin);

    // Read scattering cosines
    Nmu = static_cast<uint32_t>(ace.nxs(5) + 1);
    i = ace.jxs(5) - 1;
    for (size_t ie = 0; ie < Ne; ie++) {
      std::vector<double> cosines_for_ie = ace.xss(i, Nmu);
      i += Nmu;

      // Check cosines
      if (!std::is_sorted(cosines_for_ie.begin(), cosines_for_ie.end())) {
        std::string mssg =
            "STIncoherentElastic::STIncoherentElastic: Cosines are not sored "
            "for "
            "incoming energy index " +
            std::to_string(ie) + ".";
        throw PNDLException(mssg, __FILE__, __LINE__);
      }

      if (cosines_for_ie.front() < -1.) {
        std::string mssg =
            "STIncoherentElastic::STIncoherentElastic: Lowest cosine is less "
            "than -1 for incoming energy index " +
            std::to_string(ie) + ".";
        throw PNDLException(mssg, __FILE__, __LINE__);
      }

      if (cosines_for_ie.back() > 1.) {
        std::string mssg =
            "STIncoherentElastic::STIncoherentElastic: Largest cosine is "
            "greater "
            "than 1 for incoming eneergy index " +
            std::to_string(ie) + ".";
        throw PNDLException(mssg, __FILE__, __LINE__);
      }

      cosines_.push_back(cosines_for_ie);
    }  // For all incoming energies
  }
}

}  // namespace pndl