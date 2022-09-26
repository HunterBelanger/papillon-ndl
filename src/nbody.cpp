/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2022, Hunter Belanger
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
#include <PapillonNDL/nbody.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <cmath>

#include "constants.hpp"

namespace pndl {

NBody::NBody(const ACE& ace, std::size_t i, double iQ)
    : n_(), Ap_(), A_(), Q_(iQ) {
  n_ = ace.xss<uint32_t>(i);
  Ap_ = ace.xss(i + 1);
  A_ = ace.awr();

  if ((n_ != 3) && (n_ != 4) && (n_ != 5)) {
    std::string mssg =
        "n may only be 3, 4, or 5. Was given n = " + std::to_string(n_) +
        ". Index to XSS block is " + std::to_string(i) + ".";
    throw PNDLException(mssg);
  }

  if (Ap_ <= 0.) {
    std::string mssg =
        "Total mass ratio of all particles (Ap) must be greater than zero. "
        "Index to XSS block is " +
        std::to_string(i) + ".";
    throw PNDLException(mssg);
  }

  if (A_ <= 0.) {
    std::string mssg = "Atomic weight ratio (AWR) must be greater than zero.";
    // No need to give index here, as this value is taken directly from
    // the ACE file data.
    throw PNDLException(mssg);
  }
}

NBody::NBody(uint16_t n, double Ap, double AWR, double Q)
    : n_(n), Ap_(Ap), A_(AWR), Q_(Q) {
  if ((n_ != 3) && (n_ != 4) && (n_ != 5)) {
    std::string mssg =
        "n may only be 3, 4, or 5. Was given n = " + std::to_string(n_) + ".";
    throw PNDLException(mssg);
  }

  if (Ap_ <= 0.) {
    std::string mssg =
        "Total mass ratio of all particles (Ap) must be greater than zero.";
    throw PNDLException(mssg);
  }

  if (A_ <= 0.) {
    std::string mssg = "Atomic weight ratio (AWR) must be greater than zero.";
    throw PNDLException(mssg);
  }
}

AngleEnergyPacket NBody::sample_angle_energy(
    double E_in, std::function<double()> rng) const {
  const double Emax = ((Ap_ - 1.) / Ap_) * ((A_ / (A_ + 1.)) * E_in + Q_);
  const double x = maxwellian_spectrum(rng);
  double y = 0.;
  switch (n_) {
    case 3: {
      y = maxwellian_spectrum(rng);
      break;
    }
    case 4: {
      y = -std::log(rng() * rng() * rng());
      break;
    }
    case 5: {
      y = -std::log(rng() * rng() * rng() * rng()) -
          std::log(rng()) * std::pow(std::cos(0.5 * PI * rng()), 2.);
      break;
    }
  }

  double E_out = (x / (x + y)) * Emax;
  double mu = 2. * rng() - 1.;

  return {mu, E_out};
}

double NBody::maxwellian_spectrum(std::function<double()>& rng) const {
  const double a = PI * rng() / 2.;
  return -(std::log(rng()) + std::log(rng()) * std::cos(a) * std::cos(a));
}

std::optional<double> NBody::angle_pdf(double /*E_in*/, double /*mu*/) const {
  // NBody is isotropic in the CM frame
  return 0.5;
}

std::optional<double> NBody::pdf(double E_in, double /*mu*/,
                                 double E_out) const {
  double Emax = ((Ap_ - 1.) / Ap_) * ((A_ / (A_ + 1.)) * E_in + Q_);
  double C = 0;
  switch (n_) {
    case 3:
      C = 4. / (PI * Emax * Emax);
      break;

    case 4:
      C = 105. / (32. * std::pow(Emax, 7. / 2.));
      break;

    case 5:
      C = 256. / (14. * PI * std::pow(Emax, 5.));
      break;
  }

  double p = C * std::sqrt(E_out) *
             std::pow(Emax - E_out, (3. * static_cast<double>(n_) / 2.) - 4.);

  return p;
}

}  // namespace pndl
