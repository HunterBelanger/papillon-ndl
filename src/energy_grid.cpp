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
#include <PapillonNDL/energy_grid.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>

namespace pndl {

EnergyGrid::EnergyGrid(const ACE& ace, uint32_t NBINS)
    : energy_values_({0.}), bin_pointers_(), u_min(), du() {
  std::vector<float> raw_egrid = ace.xss<float>(ace.ESZ(), ace.nxs(2));

  energy_values_ = shared_span<float>(raw_egrid.begin(), raw_egrid.end());

  if (!std::is_sorted(energy_values_.begin(), energy_values_.end())) {
    std::string mssg =
        "EnergyGrid::EnergyGrid: Energy values are not sorted. Index in the "
        "XSS block is " +
        std::to_string(ace.ESZ()) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (energy_values_.front() <= 0.) {
    std::string mssg =
        "EnergyGrid::EnergyGrid: Nevative or zero values in energy grid. Index "
        "in the XSS block is " +
        std::to_string(ace.ESZ()) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  hash_energy_grid(NBINS);
}

EnergyGrid::EnergyGrid(const std::vector<double>& energy, uint32_t NBINS)
    : energy_values_(energy.size(), 0.), bin_pointers_(), u_min(), du() {
  if (!std::is_sorted(energy_values_.begin(), energy_values_.end())) {
    std::string mssg = "EnergyGrid::EnergyGrid: Energy values are not sorted.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (energy_values_.front() <= 0.) {
    std::string mssg =
        "EnergyGrid::EnergyGrid: Nevative or zero values in energy grid.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  for (std::size_t i = 0; i < energy.size(); i++) {
    energy_values_[i] = static_cast<float>(energy[i]);
  }

  hash_energy_grid(NBINS);
}

EnergyGrid::EnergyGrid(const std::vector<float>& energy, uint32_t NBINS)
    : energy_values_(energy.begin(), energy.end()),
      bin_pointers_(),
      u_min(),
      du() {
  if (!std::is_sorted(energy_values_.begin(), energy_values_.end())) {
    std::string mssg =
        "EnergyGrid::EnergyGrid: Energy values are not sorted.\n";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (energy_values_.front() <= 0.) {
    std::string mssg =
        "EnergyGrid::EnergyGrid: Nevative or zero values in energy grid.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  hash_energy_grid(NBINS);
}

double EnergyGrid::operator[](size_t i) const { return energy_values_[i]; }

size_t EnergyGrid::size() const { return energy_values_.size(); }

shared_span<float> EnergyGrid::grid() const { return energy_values_; }

double EnergyGrid::min_energy() const { return energy_values_.front(); }

double EnergyGrid::max_energy() const { return energy_values_.back(); }

void EnergyGrid::hash_energy_grid(uint32_t NBINS) {
  // Generate pointers for lethargy bins
  u_min = std::log(energy_values_.front());
  double u_max = std::log(energy_values_.back());
  du = (u_max - u_min) / static_cast<double>(NBINS);

  bin_pointers_.clear();
  bin_pointers_.reserve(NBINS + 1);

  double E = energy_values_.front();
  size_t i = 0;

  // Start by storing index to u_min which is 0
  bin_pointers_.push_back(0);

  // Get energy index for each lethargy bin bound
  for (size_t b = 1; b < NBINS + 1; b++) {
    E *= std::exp(du);

    i = std::distance(
            energy_values_.begin(),
            std::lower_bound(energy_values_.begin(), energy_values_.end(), E)) -
        1;

    bin_pointers_.push_back(static_cast<uint32_t>(i));
  }
}

}  // namespace pndl
