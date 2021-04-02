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
#include <PapillonNDL/region_1d.hpp>
#include <PapillonNDL/watt.hpp>
#include <cmath>

#include "constants.hpp"

namespace pndl {

Watt::Watt(const ACE& ace, size_t i) : a_(), b_(), restriction_energy_() {
  size_t original_i = i;
  uint32_t NR = ace.xss<uint32_t>(i);
  uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
  std::vector<uint32_t> NBT_a;
  std::vector<Interpolation> INT_a;

  if (NR == 0) {
    NBT_a = {NE};
    INT_a = {Interpolation::LinLin};
  } else {
    NBT_a = ace.xss<uint32_t>(i + 1, NR);
    INT_a = ace.xss<Interpolation>(i + 1 + NR, NR);
  }

  // Get energy grid
  std::vector<double> energy_a = ace.xss(i + 2 + 2 * NR, NE);
  std::vector<double> a = ace.xss(i + 2 + 2 * NR + NE, NE);

  // Create Function1D pointer
  try {
    if (NBT_a.size() == 1) {
      a_ = std::make_shared<Region1D>(energy_a, a, INT_a[0]);
    } else {
      a_ = std::make_shared<MultiRegion1D>(NBT_a, INT_a, energy_a, a);
    }
  } catch (PNDLException& error) {
    std::string mssg =
        "Watt::Watt: Could not construct Tabular1D for the 'a'. Index in the "
        "XSS block is i = " +
        std::to_string(original_i) + ".";
    error.add_to_exception(mssg, __FILE__, __LINE__);
    throw error;
  }

  // Reset i for b
  i = i + 2 + 2 * NR + 2 * NE;

  NR = ace.xss<uint32_t>(i);
  NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
  std::vector<uint32_t> NBT_b;
  std::vector<Interpolation> INT_b;

  if (NR == 0) {
    NBT_b = {NE};
    INT_b = {Interpolation::LinLin};
  } else {
    NBT_b = ace.xss<uint32_t>(i + 1, NR);
    INT_b = ace.xss<Interpolation>(i + 1 + NR, NR);
  }

  // Get energy grid
  std::vector<double> energy_b = ace.xss(i + 2 + 2 * NR, NE);
  std::vector<double> b = ace.xss(i + 2 + 2 * NR + NE, NE);

  // Create Function1D pointer
  try {
    if (NBT_b.size() == 1) {
      b_ = std::make_shared<Region1D>(energy_b, b, INT_b[0]);
    } else {
      b_ = std::make_shared<MultiRegion1D>(NBT_b, INT_b, energy_b, b);
    }
  } catch (PNDLException& error) {
    std::string mssg =
        "Watt::Watt: Could not construct Tabular1D for the 'b'. Index in the "
        "XSS block is i = " +
        std::to_string(original_i) + ".";
    error.add_to_exception(mssg, __FILE__, __LINE__);
    throw error;
  }

  // Get restriction energy
  restriction_energy_ = ace.xss(i + 2 + 2 * NR + 2 * NE);
}

Watt::Watt(std::shared_ptr<Tabulated1D> a, std::shared_ptr<Tabulated1D> b,
           double restriction_energy)
    : a_(a), b_(b), restriction_energy_(restriction_energy) {}

double Watt::sample_energy(double E_in, std::function<double()> rng) const {
  double a = (*a_)(E_in);
  double b = (*b_)(E_in);
  double w = 0.;
  double xi1 = 0.;
  double xi2 = 0.;
  double xi3 = 0.;
  double c = 0.;
  double E_out = E_in;

  bool sampled = false;
  while (!sampled) {
    xi1 = rng();
    xi2 = rng();
    xi3 = rng();

    c = std::cos(PI * xi3 / 2.);

    w = -a * (std::log(xi1) + std::log(xi2) * c * c);

    E_out = w + 0.25 * a * a * b + (2. * rng() - 1.) * std::sqrt(a * a * b * w);

    if (E_out >= 0. && E_out <= (E_in - restriction_energy_)) sampled = true;
  }

  return E_out;
}

double Watt::pdf(double E_in, double E_out) const {
  double du = E_in - restriction_energy_;
  if (E_out < 0. || E_out > du) return 0.;

  double a = (*a_)(E_in);
  double b = (*b_)(E_in);
  double I = 0.5 * std::sqrt(PI * a * a * a * b / 4.) * std::exp(a * b / 4.);
  I *= std::erf(std::sqrt(du / a) - std::sqrt(a * b / 4.)) +
       std::erf(std::sqrt(du / a) + std::sqrt(a * b / 4.));
  I -= a * std::exp(-du / a) * std::sinh(std::sqrt(b * du));
  return (std::exp(-E_out / a) / I) * std::sinh(std::sqrt(b * E_out));
}

std::shared_ptr<Tabulated1D> Watt::a() const { return a_; }

std::shared_ptr<Tabulated1D> Watt::b() const { return b_; }

double Watt::U() const { return restriction_energy_; }

}  // namespace pndl
