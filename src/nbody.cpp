/*
 * Copyright 2020, Hunter Belanger
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
#include <PapillonNDL/nbody.hpp>
#include <cmath>

#include "constants.hpp"

namespace pndl {

NBody::NBody(const ACE& ace, size_t i, double iQ) : n_(), Ap_(), A_(), Q_(iQ) {
  n_ = ace.xss<uint32_t>(i);
  Ap_ = ace.xss(i + 1);
  A_ = ace.awr();
}

AngleEnergyPacket NBody::sample_angle_energy(
    double E_in, std::function<double()> rng) const {
  double Emax = ((Ap_ + 1.) / Ap_) * ((A_ / (A_ + 1.)) * E_in + Q_);
  double x = maxwellian_spectrum(rng);
  double y = 0.;
  double xi1, xi2, xi3, xi4, xi5, xi6;
  switch (n_) {
    case 3:
      y = maxwellian_spectrum(rng);
      break;
    case 4:
      xi1 = rng();
      xi2 = rng();
      xi3 = rng();
      y = -std::log(xi1 * xi2 * xi3);
      break;
    case 5:
      xi1 = rng();
      xi2 = rng();
      xi3 = rng();
      xi4 = rng();
      xi5 = rng();
      xi6 = rng();
      y = -std::log(xi1 * xi2 * xi3 * xi4) -
          std::log(xi5) * std::pow(std::cos(0.5 * PI * xi6), 2.);
      break;
  }

  double E_out = (x / (x + y)) * Emax;
  double mu = 2. * rng() - 1.;

  return {mu, E_out};
}

double NBody::maxwellian_spectrum(std::function<double()>& rng) const {
  double xi1 = rng();
  double xi2 = rng();
  double xi3 = rng();

  double a = PI * xi3 / 2.;

  return -(std::log(xi1) + std::log(xi2) * std::cos(a) * std::cos(a));
}

uint32_t NBody::n() const { return n_; }

double NBody::Ap() const { return Ap_; }

double NBody::A() const { return A_; }

double NBody::Q() const { return Q_; }

}  // namespace pndl
