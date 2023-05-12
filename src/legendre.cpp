/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2023, Hunter Belanger
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

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/legendre.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <cmath>
#include <functional>

#include "constants.hpp"

namespace pndl {

Legendre::Legendre() : a_({0.5}), pdf_max_(0.) {
  this->find_max();

  if (this->positive() == false) {
    throw PNDLException("Legendre distribution is negative within [-1,1].");
  }
}

Legendre::Legendre(const std::vector<double>& a) : a_(a), pdf_max_(0.) {
  // Add the a_[0] = 1. moment.
  a_.insert(a_.begin(), 1.);

  // Multiply all coefficients by (2l + 1)/2.
  for (std::size_t l = 0; l < a_.size(); l++) {
    a_[l] *= (2. * static_cast<double>(l) + 1.) / 2.;
  }

  this->find_max();

  if (this->positive() == false) {
    throw PNDLException("Legendre distribution is negative within [-1,1].");
  }
}

double Legendre::pdf(double mu) const {
  double pdf = 0.;

  for (std::size_t l = 0; l < a_.size(); l++) {
    pdf += a_[l] * std::legendre(static_cast<unsigned int>(l), mu);
  }

  return pdf;
}

double Legendre::sample_mu(const std::function<double()>& rng) const {
  double mu = -2.;

  while (true) {
    mu = 2. * rng() - 1.;
    double P_accept = pdf(mu) / pdf_max_;
    if (rng() < P_accept) break;
  }

  return mu;
}

void Legendre::set_moment(std::size_t l, double a) {
  // Do not allow modifying the l=0 moment, as this must always be 0.
  if (l == 0) {
    throw PNDLException("Cannot modify legendre moment 0.");
  }

  // Make sure the legendre moment we want to access is okay.
  if (l >= a_.size()) a_.resize(l + 1, 0.);

  // Set the moment
  a_[l] = a * (2. * static_cast<double>(l) + 1.) / 2.;

  // Find new max
  this->find_max();

  // Check positivity
  if (this->positive() == false) {
    throw PNDLException("Legendre distribution is negative within [-1,1].");
  }
}

void Legendre::find_max() {
  pdf_max_ = 0.;
  double mu = -1.;
  double pdf_mu = this->pdf(mu);
  double pdf_max = pdf_mu;

  constexpr double d_mu = 0.001;
  for (std::size_t i = 0; i < 2000; i++) {
    mu += d_mu;
    pdf_mu = this->pdf(mu);
    if (pdf_mu > pdf_max) pdf_max = pdf_mu;
  }

  pdf_max_ = 1.001 * pdf_max;
}

bool Legendre::positive() const {
  constexpr double d_mu = 0.001;

  double mu = -1.;
  for (std::size_t i = 0; i <= 2000; i++) {
    if (this->pdf(mu) < 0.) return false;

    mu += d_mu;
  }

  return true;
}

}  // namespace pndl
