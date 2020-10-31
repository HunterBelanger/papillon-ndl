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
#include <PapillonNDL/energy_angle_table.hpp>
#include <PapillonNDL/interpolation.hpp>
#include <cmath>

namespace pndl {

EnergyAngleTable::EnergyAngleTable(const ACE& ace, size_t i)
    : energy_(), pdf_(), cdf_(), angles_(), interp_() {
  interp_ = ace.xss<Interpolation>(i);
  if ((interp_ != Interpolation::Histogram) &&
      (interp_ != Interpolation::LinLin)) {
    throw std::runtime_error("EnergyAngleTable: Invalid interpolation");
  }
  uint32_t NP = ace.xss<uint32_t>(i + 1);
  energy_ = ace.xss(i + 2, NP);

  pdf_ = ace.xss(i + 2 + NP, NP);
  cdf_ = ace.xss(i + 2 + NP + NP, NP);

  if (!std::is_sorted(energy_.begin(), energy_.end())) {
    throw std::runtime_error("EnergyAngleTable: Energies are not sorted");
  }

  if (!std::is_sorted(cdf_.begin(), cdf_.end())) {
    throw std::runtime_error("EnergyAngleTable: CDF is not sorted");
  }

  std::vector<int32_t> locs = ace.xss<int32_t>(i + 2 + NP + NP + NP, NP);
  for (const auto& loc : locs) {
    size_t l = ace.DLW() + std::abs(loc) - 1;
    ;
    angles_.emplace_back(ace, l);
  }
}

AngleEnergyPacket EnergyAngleTable::sample_angle_energy(
    std::function<double()> rng) const {
  double E_out, mu;
  double xi = rng();
  auto cdf_it = std::lower_bound(cdf_.begin(), cdf_.end(), xi);
  size_t l = std::distance(cdf_.begin(), cdf_it) - 1;

  if (interp_ == Interpolation::Histogram) {
    E_out = histogram_interp_energy(xi, l);
    mu = angles_[l].sample_value(rng());
    if (std::abs(mu) > 1.) mu = std::copysign(1., mu);
    return {mu, E_out};
  }

  E_out = linear_interp_energy(xi, l);

  double f = interpolation_factor(xi, cdf_[l], cdf_[l + 1]);
  if (f < 0.5)
    mu = angles_[l].sample_value(rng());
  else
    mu = angles_[l + 1].sample_value(rng());

  if (std::abs(mu) > 1.) mu = std::copysign(1., mu);

  return {mu, E_out};
}

double EnergyAngleTable::min_energy() const { return energy_.front(); }

double EnergyAngleTable::max_energy() const { return energy_.back(); }

Interpolation EnergyAngleTable::interpolation() const { return interp_; }

double EnergyAngleTable::histogram_interp_energy(double xi, size_t l) const {
  return energy_[l] + ((xi - cdf_[l]) / pdf_[l]);
}

double EnergyAngleTable::linear_interp_energy(double xi, size_t l) const {
  double m = (pdf_[l + 1] - pdf_[l]) / (energy_[l + 1] - energy_[l]);
  return energy_[l] +
         (1. / m) *
             (std::sqrt(pdf_[l] * pdf_[l] + 2. * m * (xi - cdf_[l])) - pdf_[l]);
}

const std::vector<double>& EnergyAngleTable::energy() const { return energy_; }

const std::vector<double>& EnergyAngleTable::pdf() const { return pdf_; }

const std::vector<double>& EnergyAngleTable::cdf() const { return cdf_; }

const PCTable& EnergyAngleTable::angle_table(size_t i) const {
  return angles_[i];
}

size_t EnergyAngleTable::size() const { return energy_.size(); }

}  // namespace pndl
