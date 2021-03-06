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
#include <PapillonNDL/multi_region_1d.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace pndl {

MultiRegion1D::MultiRegion1D(const std::vector<Region1D>& regions)
    : regions_(regions) {
  // Assure there are at least two regions
  if (regions_.size() < 2) {
    std::string mssg =
        "MultiRegion1D::MultiRegion1D: Must provide at least 2 regions.\n";
    mssg +=
        "Was provided with " + std::to_string(regions_.size()) + " regions.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  // Ensure all regions are ordered and not overlapping
  for (size_t i = 0; i < regions_.size(); i++) {
    if (i != regions_.size() - 1) {
      if (regions_[i].min_x() >= regions_[i + 1].min_x() ||
          regions_[i].max_x() != regions_[i + 1].min_x()) {
        // Problem with ordering
        std::string mssg =
            "MultiRegion1D::MultiRegion1D: Regions provided to MultiRegion1D "
            "constructor are\n";
        mssg += "improperly orderd.";
        throw PNDLException(mssg, __FILE__, __LINE__);
      }
    }
  }
}

MultiRegion1D::MultiRegion1D(const std::vector<uint32_t>& NBT,
                             const std::vector<Interpolation>& INT,
                             const std::vector<double>& x,
                             const std::vector<double>& y)
    : regions_() {
  // Ensure NBT and INT are the same length
  if (NBT.size() != INT.size()) {
    std::string mssg =
        "MultiRegion1D::MultiRegion1D: NBT and INT have different sizes.\n";
    mssg += "NBT.size() = " + std::to_string(NBT.size()) +
            " and INT.size() = " + std::to_string(INT.size()) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (x.size() != y.size()) {
    std::string mssg =
        "MultiRegion1D::MultiRegion1D: x and y have different sizes.\n";
    mssg += "x.size() = " + std::to_string(x.size()) +
            " and y.size() = " + std::to_string(y.size()) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  if (!std::is_sorted(x.begin(), x.end())) {
    std::string mssg = "MultiRegion1D::MultiRegion1D: x is not sorted.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  // Make 1D regions of all intervals
  size_t low = 0;
  size_t hi = 0;
  for (size_t i = 0; i < NBT.size(); i++) {
    hi = NBT[i];

    try {
      regions_.push_back(Region1D({x.begin() + low, x.begin() + hi},
                                  {y.begin() + low, y.begin() + hi}, INT[i]));
    } catch (PNDLException& error) {
      std::string mssg =
          "MultiRegion1D::MultiRegion1D: The i = " + std::to_string(i) +
          " Region1D could not be constructed\n";
      mssg += "when building MultiRegion1D.";
      error.add_to_exception(mssg, __FILE__, __LINE__);
      throw error;
    }

    low = hi - 1;

    // Check for discontinuity at region boundary
    if (low < x.size() - 1) {
      if (x[low] == x[low + 1]) low++;
    }
  }
}

double MultiRegion1D::operator()(double x) const {
  if (x <= min_x())
    return (regions_.front())(x);
  else if (x >= max_x())
    return (regions_.back())(x);

  // Get region which contains x
  auto region_it = std::lower_bound(regions_.begin(), regions_.end(), x);
  region_it--;

  return (*region_it)(x);
}

double MultiRegion1D::integrate(double x_low, double x_hi) const {
  bool inverted = x_low > x_hi;
  if (inverted) {
    double x_low_tmp = x_low;
    x_low = x_hi;
    x_hi = x_low_tmp;
  }

  // Integration may only be carried out over the function's valid domain
  if (x_low <= min_x())
    x_low = min_x();
  else if (x_low >= max_x())
    x_low = max_x();

  if (x_hi >= max_x())
    x_hi = max_x();
  else if (x_hi <= min_x())
    x_hi = min_x();

  // Get region which contains x_low
  auto region = std::lower_bound(regions_.begin(), regions_.end(), x_low);
  if (region->min_x() > x_low) region--;

  double integral = 0.;
  double x_low_lim = x_low;
  double x_upp_lim = x_hi;
  bool integrating = true;
  while (integrating) {
    if (x_low_lim < region->min_x()) x_low_lim = region->min_x();
    if (x_upp_lim > region->max_x()) x_upp_lim = region->max_x();

    integral += region->integrate(x_low_lim, x_upp_lim);

    if (x_upp_lim == x_hi)
      integrating = false;
    else {
      x_low_lim = x_upp_lim;
      x_upp_lim = x_hi;
      region++;
    }
  }

  if (inverted) integral *= -1.;

  return integral;
}

const Region1D& MultiRegion1D::operator[](size_t i) const {
  return regions_[i];
}

std::vector<uint32_t> MultiRegion1D::breakpoints() const {
  std::vector<uint32_t> brks;

  for (const auto& r : regions_) {
    uint32_t nr = static_cast<uint32_t>(r.breakpoints()[0]);
    if (!brks.empty()) nr += brks.back();
    brks.push_back(nr);
  }

  return brks;
}

std::vector<Interpolation> MultiRegion1D::interpolation() const {
  std::vector<Interpolation> interps;

  for (const auto& r : regions_) {
    interps.push_back(r.interpolation()[0]);
  }

  return interps;
}

std::vector<double> MultiRegion1D::x() const {
  std::vector<double> x_;

  for (const auto& r : regions_) {
    std::vector<double> x_r = r.x();
    x_.insert(x_.end(), x_r.begin(), x_r.end());
  }

  return x_;
}

std::vector<double> MultiRegion1D::y() const {
  std::vector<double> y_;

  for (const auto& r : regions_) {
    std::vector<double> y_r = r.y();
    y_.insert(y_.end(), y_r.begin(), y_r.end());
  }

  return y_;
}

size_t MultiRegion1D::size() const { return regions_.size(); }

double MultiRegion1D::min_x() const { return regions_.front().min_x(); }

double MultiRegion1D::max_x() const { return regions_.back().max_x(); }

}  // namespace pndl
