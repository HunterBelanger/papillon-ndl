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
#include <PapillonNDL/st_coherent_elastic.hpp>

namespace pndl {

STCoherentElastic::STCoherentElastic(const ACE& ace)
    : bragg_edges_(), structure_factor_sum_() {
  // Fist make sure ACE file does indeed give coherent elastic scattering
  int32_t elastic_mode = ace.nxs(4);
  if (elastic_mode != 4) {
    std::string mssg =
        "STCoherentElastic::STCoherentElastic: Provided ACE file does not have "
        "coherent elastic scattering.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  // Get index to Bragg edge and structure data
  int32_t i = ace.jxs(3) - 1;
  uint32_t Ne = ace.xss<uint32_t>(i);
  bragg_edges_ = ace.xss(i + 1, Ne);
  structure_factor_sum_ = ace.xss(i + 1 + Ne, Ne);

  // Make sure Bragg edges are all positive and sorted
  if (!std::is_sorted(bragg_edges_.begin(), bragg_edges_.end())) {
    std::string mssg =
        "STCoherentElastic::STCoherentElastic: Bragg edges are not sorted.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (bragg_edges_.front() < 0.) {
    std::string mssg =
        "STCoherentElastic::STCoherentElastic: Negative Bragg edges found.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }
}

}  // namespace pndl