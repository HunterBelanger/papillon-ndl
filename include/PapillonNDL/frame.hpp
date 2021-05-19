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
#ifndef PAPILLON_NDL_FRAME_H
#define PAPILLON_NDL_FRAME_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/angle_energy.hpp>
#include <cmath>
#include <cstdint>

namespace pndl {

/**
 * @brief Enum to indicate the frame of reference for secondary data
 *        angle and energy data.
 */
enum class Frame : uint32_t {
  Lab = 1, /**< Laboratory frame */
  CM = 2,  /**< Center of Mass frame */
};

/**
 * @brief A struct contianing helper methods to convert scattering angle and
 *        energies provided in the center of mass frame, to the lab frame.
 */
struct CMToLab {
  /**
   * @brief Transfrom mu and Eout from the CM frame to the Lab frame.
   * @param Ein Incident energy of the particle.
   * @param A Atomic weight ratio of the target nuclide.
   * @param mu Scattering angle in the center of mass frame. The value
   *           is changed to the scattering angle in the lab frame
   *           upon return.
   * @param Eout Scattering energy in the center of mass frame. The value
   *             is changed to the scattering energy in the lab frame
   *             upon return.
   */
  static void transform(double Ein, double A, double& mu, double& Eout) {
    double Eout_lab =
        Eout + (Ein + 2. * mu * (A + 1.) * std::sqrt(Ein * Eout)) /
                   std::pow(A + 1., 2.);

    mu = mu * std::sqrt(Eout / Eout_lab) +
         (1. / (A + 1.)) * std::sqrt(Ein / Eout_lab);

    Eout = Eout_lab;
  }

  /**
   * @brief Transfrom an AngleEnergyPacket from the CM frame to the Lab frame.
   * @param Ein Incident energy of the particle.
   * @param A Atomic weight ratio of the target nuclide.
   * @param ae AngleEnergyPacket which initially contains the scattering angle
   *           and energy in the CM frame, which is changed to the lab frame
   *           upon return.
   */
  static void transform(double Ein, double A, AngleEnergyPacket& ae) {
    double E_cm = ae.energy;
    double mu_cm = ae.cosine_angle;

    double E_lab =
        E_cm + (Ein + 2. * mu_cm * (A + 1.) * std::sqrt(Ein * E_cm)) /
                   std::pow(A + 1., 2.);

    double mu_lab = mu_cm * std::sqrt(E_cm / E_lab) +
                    (1. / (A + 1.)) * std::sqrt(Ein / E_lab);

    ae.cosine_angle = mu_lab;
    ae.energy = E_lab;
  }
};

}  // namespace pndl

#endif
