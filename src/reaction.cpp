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
#include <PapillonNDL/reaction.hpp>
#include <memory>

/**
 * @file
 * @author Hunter Belanger
 */

namespace pndl {

Reaction<CrossSection>::Reaction(const ACE& ace, std::size_t indx,
                                 std::shared_ptr<EnergyGrid> egrid)
    : ReactionBase(ace, indx), xs_(nullptr) {
  try {
    uint32_t loca = ace.xss<uint32_t>(ace.LSIG() + indx);
    xs_ = std::make_shared<CrossSection>(ace, ace.SIG() + loca - 1, egrid);
    threshold_ = xs_->energy(0);
  } catch (PNDLException& error) {
    std::string mssg = "Could not create cross section for MT = " +
                       std::to_string(this->mt()) + ".";
    error.add_to_exception(mssg);
    throw error;
  }
}

Reaction<CrossSection>::Reaction(const ACE& ace, std::size_t indx,
                                 std::shared_ptr<EnergyGrid> egrid,
                                 const Reaction& reac)
    : ReactionBase(reac), xs_(nullptr) {
  // make sure the MT values agree
  if (this->mt() != reac.mt()) {
    std::string mssg = "MT from ACE file doesn't match MT from reaction.";
    throw PNDLException(mssg);
  }

  // Get XS from new ACE
  try {
    uint32_t loca = ace.xss<uint32_t>(ace.LSIG() + indx);
    xs_ = std::make_shared<CrossSection>(ace, ace.SIG() + loca - 1, egrid);
    threshold_ = xs_->energy(0);
  } catch (PNDLException& error) {
    std::string mssg =
        "Could not create cross section for MT = " + std::to_string(mt_) + ".";
    error.add_to_exception(mssg);
    throw error;
  }
}

}  // namespace pndl
