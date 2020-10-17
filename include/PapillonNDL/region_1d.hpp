/*
 * Copyright 2020, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * Ce  logiciel  est  un  programme  informatique  servant  à  résoudre
 * l'équation de Boltzmann  pour les  neutrons, avec  des énergies continus,
 * en trois  dimensions  spatiales, utilisant la méthode  Monte Carlo.
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
#ifndef PAPILLON_NDL_REGION_1D_H
#define PAPILLON_NDL_REGION_1D_H

#include <PapillonNDL/interpolation.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <vector>

namespace pndl {

class Region1D : public Tabulated1D {
 public:
  Region1D(const std::vector<double>& i_x, const std::vector<double>& i_y,
           Interpolation interp = Interpolation::LinLin);
  ~Region1D() = default;

  // Methods required by Function1D
  double operator()(double x) const override final;
  double integrate(double x_low, double x_hi) const override final;

  // Methods required by Tabulated1D
  std::vector<uint32_t> breakpoints() const override final;
  std::vector<Interpolation> interpolation() const override final;
  std::vector<double> x() const override final;
  std::vector<double> y() const override final;

  // Extra Methods
  size_t size() const;
  double min_x() const;
  double max_x() const;

 private:
  Interpolation interpolation_;
  std::vector<double> x_;
  std::vector<double> y_;
};

// This operator overload is provided only to accomodate the
// std::lower_bound algorithm, in the MultiRegion1D class
bool operator<(const Region1D& R, const double& X);

}  // namespace pndl

#endif
