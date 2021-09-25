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
#ifndef PAPILLON_NDL_FRAME_H
#define PAPILLON_NDL_FRAME_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/angle_energy.hpp>
#include <cmath>
#include <cstdint>
#include <optional>

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
    CMToLab::transform(Ein, A, ae.cosine_angle, ae.energy);
  }

  /**
   * @brief Computes the Jacobian to transform the angular PDF from the center
   *        of mass frame, to the lab frame. This quantity is dmu_cm/dmu_lab.
   *        This formula comes from the public MCNP theory manual.
   * @param Ein Incident energy in the lab frame.
   * @param A Atomic weight ratio of the target nuclide.
   * @param mu Cosine of the scattering angle in the lab frame.
   * @param Eout Scattering energy in the lab frame.
   */
  static double angle_jacobian(double Ein, double A, double mu, double Eout) {
    double C = Ein + 2. * (A + 1.) * std::sqrt(Ein * Eout) *
                         (mu - (1. / (A + 1.)) * std::sqrt(Ein / Eout));
    double Eout_cm = Eout - (C / ((A + 1.) * (A + 1.)));
    return std::sqrt(Eout / Eout_cm) /
           (1. - (mu / (A + 1.)) * std::sqrt(Ein / Eout));
  }

  /**
   * @brief Computes the Jacobian to transform the angular PDF from the center
   *        of mass frame, to the lab frame. This quantity is dmu_cm/dmu_lab.
   *        This formula comes from the public MCNP theory manual.
   * @param Ein Incident energy in the lab frame.
   * @param A Atomic weight ratio of the target nuclide.
   * @param ae AngleEnergyPacket which contains the scattering angle and energy
   *           in the lab frame.
   */
  static double angle_jacobian(double Ein, double A, AngleEnergyPacket ae) {
    double mu = ae.cosine_angle;
    double Eout = ae.energy;
    double C = Ein + 2. * (A + 1.) * std::sqrt(Ein * Eout) *
                         (mu - (1. / (A + 1.)) * std::sqrt(Ein / Eout));
    double Eout_cm = Eout - (C / ((A + 1.) * (A + 1.)));
    return std::sqrt(Eout / Eout_cm) /
           (1. - (mu / (A + 1.)) * std::sqrt(Ein / Eout));
  }

  /**
   * @brief Computes the Jacobian to transform the angular PDF from the center
   *        of mass frame, to the lab frame. This quantity is dmu_cm/dmu_lab.
   * @param Ein Incident energy in the lab frame.
   * @param A Atomic weight ratio of the target nuclide.
   * @param Q Q-value of the reaction.
   * @param mu Cosine of the scattering angle in the lab frame.
   * @param mu_cm Cosine of the scattering angle in the center of mass frame.
   */
  static double angle_jacobian(double Ein, double A, double Q, double mu,
                               double mu_cm) {
    // This form is derived from [1].
    double Ec = A * Ein / (A + 1.);
    double g = std::sqrt((1. / (A * A)) * (Ec / (Ec + Q)));

    double numerator =
        std::pow(g + mu_cm, 2.) * std::sqrt(1. - (mu_cm * mu_cm));
    double denominator = mu * mu * (1. + g * mu_cm) * std::sqrt(1. - mu * mu);

    return numerator / denominator;
  }

  /**
   * @brief Computes the Jacobian to transform the joint PDF from the center
   *        of mass frame, to the lab frame. This quantity is
   *        (dmu_cm/dmu_lab) * (dE_cm/dE_lab).
   * @param Eout Scattering energy in the lab frame.
   * @param Eout_cm Scattering energy in the center of mass frame.
   */
  static double jacobian(double Eout, double Eout_cm) {
    return std::sqrt(Eout / Eout_cm);
  }
};

/**
 * @brief A struct contianing helper methods to convert scattering angle and
 *        energies provided in the lab frame, to the center of mass frame.
 */
struct LabToCM {
  /**
   * @brief Transfrom mu and Eout from the Lab frame to the CM frame.
   * @param Ein Incident energy of the particle.
   * @param A Atomic weight ratio of the target nuclide.
   * @param mu Scattering angle in the lab frame. The value is changed to the
   *           scattering angle in the center of mass frame upon return.
   * @param Eout Scattering energy in the lab frame. The value is changed to
   *             the scattering energy in the center of mass frame upon return.
   */
  static void transform(double Ein, double A, double& mu, double& Eout) {
    double C = Ein + 2. * (A + 1.) * std::sqrt(Ein * Eout) *
                         (mu - (1. / (A + 1.)) * std::sqrt(Ein / Eout));
    double Eout_cm = Eout - (C / ((A + 1.) * (A + 1.)));
    double mu_cm = std::sqrt(Eout / Eout_cm) *
                   (mu - (1. / (A + 1.)) * std::sqrt(Ein / Eout));

    mu = mu_cm;
    Eout = Eout_cm;
  }

  /**
   * @brief Transfrom an AngleEnergyPacket from the Lab frame to the CM frame.
   * @param Ein Incident energy of the particle.
   * @param A Atomic weight ratio of the target nuclide.
   * @param ae AngleEnergyPacket which initially contains the scattering angle
   *           and energy in the Lab frame, which is changed to the center of
   *           mass frame upon return.
   */
  static void transform(double Ein, double A, AngleEnergyPacket& ae) {
    LabToCM::transform(Ein, A, ae.cosine_angle, ae.energy);
  }

  /**
   * @brief Calculates all possible values for the scattering angle in the
   *        center of mass frame, given a scattering angle in the lab frame.
   * @param Ein Incident energy of the particle.
   * @param A Atomic weight ratio of the target nuclide.
   * @param Q Q-value of the reaction.
   * @param mu Cosine of the scattering angle in the lab frame.
   */
  static std::pair<std::optional<double>, std::optional<double>> angle(
      double Ein, double A, double Q, double mu) {
    // This version does not require the outgoing energy in the lab frame, as
    // it uses the relation provided by Eq. 2 and Eq. 3 in [1].

    double Ec = A * Ein / (A + 1.);
    double g = std::sqrt((1. / (A * A)) * (Ec / (Ec + Q)));
    double mu_sqr = mu * mu;
    double n = (1. - mu_sqr) / mu_sqr;
    double a = 1. + n;
    double b = 2. * g * n;
    double c = n * g * g - 1.;

    // We have an equation for mu_cm which takes the form
    // a * mu_cm^2 + b * mu_cm + c = 0
    // This equation can sometimes have two solutions. According to [2],
    // we should use the largest solution always ?

    double b2_4ac = b * b - 4. * a * c;

    if (b2_4ac < 0.) {
      // Not possible to scatter with angle mu in the lab frame.
      return {std::nullopt, std::nullopt};
    } else if (mu == 0.) {
      // A singularity exists at mu = 0.. We treat this case first.
      // The center of the parabola exists at mu_cm = -gn/(1+n). In the limit
      // of u -> 0, we have the limit n -> inf, which produces mu_cm = -g
      return {-g, std::nullopt};
    } else if (b2_4ac == 0.) {
      // One possible solution
      double a1 = -2. * c / b;

      if (-1. <= a1 || a1 <= 1.)
        return {a1, std::nullopt};
      else
        return {std::nullopt, std::nullopt};
    } else {
      // Two possible solutions. Both SHOULD be valid, but we check anyway
      double a1 = (-b + std::sqrt(b2_4ac)) / (2. * a);
      double a2 = 2. * c / (-b + std::sqrt(b2_4ac));

      std::optional<double> opt1 = std::nullopt;
      std::optional<double> opt2 = std::nullopt;

      if (-1. <= a1 || a1 <= 1.) opt1 = a1;
      if (-1. <= a2 || a2 <= 1.) opt2 = a2;

      return {opt1, opt2};
    }
  }
};

}  // namespace pndl

// References
//
//  [1] P. F. Zweifel and H. Hurwitz, “Transformation of Scattering Cross
//      Sections,” J Appl Phys, vol. 25, no. 10, pp. 1241–1245, 1954,
//      doi: 10.1063/1.1721536.
//
//  [2] KINEMATICS II: A NONRELATIVISTIC KINEMATICS FORTRAN PROGRAM TO AID
//      ANALYSIS OF NUCLEAR REACTION ANGULAR DISTRIBUTION DATA, ORNL-3251

#endif
