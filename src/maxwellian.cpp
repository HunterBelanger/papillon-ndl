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
#include <PapillonNDL/maxwellian.hpp>
#include <PapillonNDL/multi_region_1d.hpp>
#include <PapillonNDL/region_1d.hpp>
#include <cmath>

#include "constants.hpp"

namespace pndl {

Maxwellian::Maxwellian(const ACE& ace, size_t i)
    : temperature_(), restriction_energy_() {
  uint32_t NR = ace.xss<uint32_t>(i);
  uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
  std::vector<uint32_t> NBT;
  std::vector<Interpolation> INT;

  if (NR == 0) {
    NBT = {NE};
    INT = {Interpolation::LinLin};
  } else {
    NBT = ace.xss<uint32_t>(i + 1, NR);
    INT = ace.xss<Interpolation>(i + 1 + NR, NR);
  }

  // Get energy grid, convert from MeV to eV
  std::vector<double> energy = ace.xss(i + 2 + 2 * NR, NE);
  for (auto& E : energy) E *= MEV_TO_EV;

  std::vector<double> temperature = ace.xss(i + 2 + 2 * NR + NE, NE);
  // Change temperatures from MeV to eV
  for (auto& T : temperature) T *= MEV_TO_EV;

  // Get restriction energy
  restriction_energy_ = ace.xss(i + 2 + 2 * NR + 2 * NE);

  // Create Function1D pointer
  if (NBT.size() == 1) {
    temperature_ = std::make_unique<Region1D>(energy, temperature, INT[0]);
  } else {
    temperature_ =
        std::make_unique<MultiRegion1D>(NBT, INT, energy, temperature);
  }
}

double Maxwellian::sample_energy(double E_in,
                                 std::function<double()> rng) const {
  double T = (*temperature_)(E_in);
  double xi1 = 0.;
  double xi2 = 0.;
  double xi3 = 0.;
  double a = 0.;
  double E_out = E_in;
  bool sampled = false;
  while (!sampled) {
    xi1 = rng();
    xi2 = rng();
    xi3 = rng();

    a = PI * xi3 / 2.;

    E_out = -T * (std::log(xi1) + std::log(xi2) * std::cos(a) * std::cos(a));

    if (E_out >= 0. && E_out <= (E_in - restriction_energy_)) sampled = true;
  }

  return E_out;
}

}  // namespace pndl
