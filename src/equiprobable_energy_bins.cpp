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
#include <PapillonNDL/equiprobable_energy_bins.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <cmath>

namespace pndl {

EquiprobableEnergyBins::EquiprobableEnergyBins(const ACE& ace, size_t i)
    : incoming_energy_(), bin_sets_() {
  // Get number of interpolation points
  uint32_t NR = ace.xss<uint32_t>(i);
  // Get number of energy points
  uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);

  // Breakpoints and interpolations are not read, as linear-linear
  // interpolation is always used between incoming energies.

  // Read energies
  incoming_energy_ = ace.xss(i + 2 + 2 * NR, NE);

  if (!std::is_sorted(incoming_energy_.begin(), incoming_energy_.end())) {
    std::string mssg = "EquiprobableEnergyBins::EquiprobableEnergyBins: ";
    mssg += "Incoming energy grid is not sorted.\n";
    mssg += "Index of EquiprobableEnergyBins in XSS block is ";
    mssg += std::to_string(i) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  uint32_t NET = ace.xss<uint32_t>(i + 2 + 2 * NR + NE);

  // Read energy bins
  for (size_t j = 0; j < NE; j++) {
    bin_sets_.push_back(ace.xss(i + 3 + 2 * NR + NE + j * NET, NET));
  }

  // Make sure that each bin set is sorted
  for (size_t j = 0; j < bin_sets_.size(); j++) {
    if (!std::is_sorted(bin_sets_[j].begin(), bin_sets_[j].end())) {
      std::string mssg = "EquiprobableEnergyBins::EquiprobableEnergyBins: ";
      mssg += std::to_string(j) + "th bin bounds are not sorted.\n";
      mssg += "Index of EquiprobableEnergyBins in XSS block is ";
      mssg += std::to_string(i) + ".";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }
}

EquiprobableEnergyBins::EquiprobableEnergyBins(
    const std::vector<double>& incoming_energy,
    const std::vector<std::vector<double>>& bin_bounds)
    : incoming_energy_(incoming_energy), bin_sets_(bin_bounds) {
  if (!std::is_sorted(incoming_energy_.begin(), incoming_energy_.end())) {
    std::string mssg = "EquiprobableEnergyBins::EquiprobableEnergyBins: ";
    mssg += "Incoming energy grid is not sorted.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (incoming_energy_.size() != bin_sets_.size()) {
    std::string mssg = "EquiprobableEnergyBins::EquiprobableEnergyBins: ";
    mssg += "The number of\nincoming energies does not match the numer ";
    mssg += "of bin sets provided.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  // Make sure that each bin set is sorted
  for (size_t j = 0; j < bin_sets_.size(); j++) {
    if (!std::is_sorted(bin_sets_[j].begin(), bin_sets_[j].end())) {
      std::string mssg = "EquiprobableEnergyBins::EquiprobableEnergyBins: ";
      mssg += std::to_string(j) + "th bin bounds are not sorted.";
      throw PNDLException(mssg, __FILE__, __LINE__);
    }
  }
}

double EquiprobableEnergyBins::sample_energy(
    double E_in, std::function<double()> rng) const {
  // Determine the index of the bounding tabulated incoming energies
  auto in_E_it =
      std::lower_bound(incoming_energy_.begin(), incoming_energy_.end(), E_in);
  if (in_E_it == incoming_energy_.begin()) {
    return sample_bins(rng(), rng(), bin_sets_.front());
  } else if (in_E_it == incoming_energy_.end()) {
    return sample_bins(rng(), rng(), bin_sets_.back());
  }

  size_t l = std::distance(incoming_energy_.begin(), in_E_it);
  l--;

  double f = (E_in - incoming_energy_[l]) /
             (incoming_energy_[l + 1] - incoming_energy_[l]);

  if (rng() > f) {
    return sample_bins(rng(), rng(), bin_sets_[l]);
  } else {
    return sample_bins(rng(), rng(), bin_sets_[l + 1]);
  }
}

double EquiprobableEnergyBins::sample_bins(
    double xi1, double xi2, const std::vector<double>& bounds) const {
  size_t bin = static_cast<size_t>(std::floor(bounds.size() * xi1));
  return (bounds[bin + 1] - bounds[bin]) * xi2 + bounds[bin];
}

const std::vector<double>& EquiprobableEnergyBins::incoming_energy() const {
  return incoming_energy_;
}

const std::vector<double>& EquiprobableEnergyBins::bin_bounds(size_t i) const {
  return bin_sets_[i];
}

size_t EquiprobableEnergyBins::size() const { return incoming_energy_.size(); }

}  // namespace pndl
