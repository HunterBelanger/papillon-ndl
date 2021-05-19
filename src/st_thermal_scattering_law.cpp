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
#include <PapillonNDL/st_thermal_scattering_law.hpp>

namespace pndl {

STThermalScatteringLaw::STThermalScatteringLaw(const ACE& ace,
                                               bool unit_based_interpolation)
    : zaid_(ace.zaid()),
      awr_(ace.awr()),
      temperature_(ace.temperature()),
      has_coherent_elastic_(false),
      has_incoherent_elastic_(false),
      coherent_elastic_(nullptr),
      incoherent_elastic_(nullptr),
      incoherent_inelastic_(nullptr) {
  // First try to read Incoherent Inelastic data, as all thermal scattering
  // laws have this
  try {
    incoherent_inelastic_ =
        std::make_shared<STIncoherentInelastic>(ace, unit_based_interpolation);
  } catch (PNDLException& err) {
    std::string mssg =
        "STThermalScatteringLaw::STThermalScatteringLaw: Could not construct "
        "Incoherent Inelastic scattering data.";
    err.add_to_exception(mssg, __FILE__, __LINE__);
    throw err;
  }

  // Currently, there MAY be EITHER Coherent or Incoherent Elastic
  // scattering. We need to check which and get the right one.
  if (ace.jxs(3) != 0) {
    // Check if we have coherent or incoherent
    int32_t elastic_mode = ace.nxs(4);
    if (elastic_mode == 4) {
      // Coherent
      try {
        coherent_elastic_ = std::make_shared<STCoherentElastic>(ace);
      } catch (PNDLException& err) {
        std::string mssg =
            "STThermalScatteringLaw::STThermalScatteringLaw: Could not "
            "construct Coherent Elastic scattering data.";
        err.add_to_exception(mssg, __FILE__, __LINE__);
        throw err;
      }
      has_coherent_elastic_ = true;
      incoherent_elastic_ = std::make_shared<STIncoherentElastic>(ace);
    } else {
      // Incoherent
      try {
        incoherent_elastic_ = std::make_shared<STIncoherentElastic>(ace);
      } catch (PNDLException& err) {
        std::string mssg =
            "STThermalScatteringLaw::STThermalScatteringLaw: Could not "
            "construct Incoherent Elastic scattering data.";
        err.add_to_exception(mssg, __FILE__, __LINE__);
        throw err;
      }
      has_incoherent_elastic_ = true;
      coherent_elastic_ = std::make_shared<STCoherentElastic>(ace);
    }
  } else {
    incoherent_elastic_ = std::make_shared<STIncoherentElastic>(ace);
    coherent_elastic_ = std::make_shared<STCoherentElastic>(ace);
  }
}
}  // namespace pndl
