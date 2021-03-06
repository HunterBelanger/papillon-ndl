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
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/region_1d.hpp>
#include <algorithm>

namespace pndl {

Region1D::Region1D(const std::vector<double>& i_x,
                   const std::vector<double>& i_y, Interpolation interp)
    : x_(i_x), y_(i_y), interpolation_(interp), interpolator() {
  if (x_.size() != y_.size()) {
    std::string mssg =
        "Region1D::Region1D: x and y have different sizes. x.size() = " +
        std::to_string(x_.size()) +
        " and y.size() = " + std::to_string(y_.size()) + ".";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  // Ensure x_ is ordered
  if (!std::is_sorted(x_.begin(), x_.end())) {
    std::string mssg = "Region1D::Region1D: x is not sorted.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  // Set interpolator
  if (interpolation_ == Interpolation::Histogram) {
    interpolator = Histogram();
  } else if (interpolation_ == Interpolation::LinLin) {
    interpolator = LinLin();
  } else if (interpolation_ == Interpolation::LinLog) {
    interpolator = LinLog();
  } else if (interpolation_ == Interpolation::LogLin) {
    interpolator = LogLin();
  } else if (interpolation_ == Interpolation::LogLog) {
    interpolator = LogLog();
  }

  auto verify_x_grid = [&i_x](auto& interp) {
    return interp.verify_x_grid(i_x.begin(), i_x.end());
  };
  auto verify_y_grid = [&i_y](auto& interp) {
    return interp.verify_y_grid(i_y.begin(), i_y.end());
  };

  std::visit(verify_x_grid, interpolator);
  std::visit(verify_y_grid, interpolator);
}

bool operator<(const Region1D& R, const double& X) { return R.min_x() < X; }

}  // namespace pndl
