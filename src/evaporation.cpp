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
#include <PapillonNDL/evaporation.hpp>
#include <PapillonNDL/multi_region_1d.hpp>
#include <PapillonNDL/region_1d.hpp>
#include <cmath>

namespace pndl {

Evaporation::Evaporation(const ACE& ace, std::size_t i)
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

  // Get energy grid
  std::vector<double> energy = ace.xss(i + 2 + 2 * NR, NE);

  std::vector<double> temperature = ace.xss(i + 2 + 2 * NR + NE, NE);

  // Get restriction energy
  restriction_energy_ = ace.xss(i + 2 + 2 * NR + 2 * NE);

  // Create Function1D pointer
  try {
    if (NBT.size() == 1) {
      temperature_ = std::make_shared<Region1D>(energy, temperature, INT[0]);
    } else {
      temperature_ =
          std::make_shared<MultiRegion1D>(NBT, INT, energy, temperature);
    }
  } catch (PNDLException& error) {
    std::string mssg =
        "Evaporation::Evaporation: Could not construct Tabular1D for the "
        "effective nuclear temperature. Index in the XSS block is i = " +
        std::to_string(i) + ".";
    error.add_to_exception(mssg, __FILE__, __LINE__);
    throw error;
  }
}

Evaporation::Evaporation(std::shared_ptr<Tabulated1D> temperature,
                         double restriction_energy)
    : temperature_(temperature), restriction_energy_(restriction_energy) {}

double Evaporation::sample_energy(double E_in,
                                  std::function<double()> rng) const {
  double T = (*temperature_)(E_in);
  double xi1 = 0.;
  double xi2 = 0.;
  double g = 0.;
  double w = 0.;
  double E_out = E_in;
  bool sampled = false;
  while (!sampled) {
    xi1 = rng();
    xi2 = rng();

    w = (E_in - restriction_energy_) / T;
    g = 1. - std::exp(-w);

    E_out = -T * std::log((1. - g * xi1) * (1. - g * xi2));

    if (E_out >= 0. && E_out <= (E_in - restriction_energy_)) sampled = true;
  }

  return E_out;
}

double Evaporation::pdf(double E_in, double E_out) const {
  double du = E_in - restriction_energy_;
  if (E_out < 0. || E_out > du) return 0.;

  double T = (*temperature_)(E_in);
  double I = T * T * (1. - std::exp(-du / T) * (1. + (du / T)));
  return (E_out / I) * std::exp(-E_out / T);
}

}  // namespace pndl
