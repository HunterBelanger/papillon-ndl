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
#include <pybind11/pybind11.h>

#include <PapillonNDL/version.hpp>

namespace py = pybind11;

extern void init_ACE(py::module&);
extern void init_XSPacket(py::module&);
extern void init_Interpolation(py::module&);
extern void init_Interpolator(py::module&);
extern void init_Frame(py::module&);
extern void init_Function1D(py::module&);
extern void init_Constant(py::module&);
extern void init_Polynomial1D(py::module&);
extern void init_Tabulated1D(py::module&);
extern void init_Sum1D(py::module&);
extern void init_Difference1D(py::module&);
extern void init_Linearize(py::module&);
extern void init_AngleLaw(py::module&);
extern void init_Isotropic(py::module&);
extern void init_EquiprobableAngleBins(py::module&);
extern void init_AngleTable(py::module&);
extern void init_Legendre(py::module&);
extern void init_AngleDistribution(py::module&);
extern void init_EnergyLaw(py::module&);
extern void init_DiscretePhoton(py::module&);
extern void init_TabularEnergy(py::module&);
extern void init_GeneralEvaporation(py::module&);
extern void init_EquiprobableEnergyBins(py::module&);
extern void init_LevelInelasticScatter(py::module&);
extern void init_Evaporation(py::module&);
extern void init_Maxwellian(py::module&);
extern void init_Watt(py::module&);
extern void init_AngleEnergyPacket(py::module&);
extern void init_AngleEnergy(py::module&);
extern void init_Uncorrelated(py::module&);
extern void init_Elastic(py::module&);
extern void init_NBody(py::module&);
extern void init_KalbachTable(py::module&);
extern void init_Kalbach(py::module&);
extern void init_DiscreteCosinesEnergies(py::module&);
extern void init_ContinuousEnergyDiscreteCosines(py::module& m);
extern void init_DirectSab(py::module& m);
extern void init_MultipleDistribution(py::module& m);
extern void init_SummedFissionSpectrum(py::module& m);
extern void init_CMDistribution(py::module& m);
extern void init_Absorption(py::module& m);
extern void init_STTSLReaction(py::module& m);
extern void init_STIncoherentInelastic(py::module& m);
extern void init_STCoherentElastic(py::module& m);
extern void init_STInoherentElasticACE(py::module& m);
extern void init_STInoherentElastic(py::module& m);
extern void init_STThermalScatteringLaw(py::module& m);
extern void init_URRPtable(py::module& m);
extern void init_PCTable(py::module&);
extern void init_EnergyAngleTable(py::module&);
extern void init_TabularEnergyAngle(py::module&);
extern void init_EnergyGrid(py::module&);
extern void init_CrossSection(py::module&);
extern void init_DelayedFamily(py::module&);
extern void init_Fission(py::module&);
extern void init_ReactionBase(py::module& m);
extern void init_STReaction(py::module& m);
extern void init_STNeutron(py::module& m);
extern void init_PRNG(py::module&);
extern void init_ZAID(py::module&);
extern void init_Element(py::module&);
extern void init_Isotope(py::module&);
extern void init_Nuclide(py::module&);
extern void init_NDLibrary(py::module&);
extern void init_MCNPLibrary(py::module&);
extern void init_SerpentLibrary(py::module&);

PYBIND11_MODULE(pyPapillonNDL, m) {
  init_ACE(m);
  init_XSPacket(m);
  init_Interpolation(m);
  init_Interpolator(m);
  init_Frame(m);
  init_Function1D(m);
  init_Constant(m);
  init_Polynomial1D(m);
  init_Tabulated1D(m);
  init_Sum1D(m);
  init_Difference1D(m);
  init_Linearize(m);
  init_AngleLaw(m);
  init_Isotropic(m);
  init_EquiprobableAngleBins(m);
  init_AngleTable(m);
  init_Legendre(m);
  init_AngleDistribution(m);
  init_EnergyLaw(m);
  init_DiscretePhoton(m);
  init_TabularEnergy(m);
  init_GeneralEvaporation(m);
  init_EquiprobableEnergyBins(m);
  init_LevelInelasticScatter(m);
  init_Evaporation(m);
  init_Maxwellian(m);
  init_Watt(m);
  init_AngleEnergyPacket(m);
  init_AngleEnergy(m);
  init_Uncorrelated(m);
  init_Elastic(m);
  init_NBody(m);
  init_KalbachTable(m);
  init_Kalbach(m);
  init_DiscreteCosinesEnergies(m);
  init_ContinuousEnergyDiscreteCosines(m);
  init_DirectSab(m);
  init_MultipleDistribution(m);
  init_SummedFissionSpectrum(m);
  init_CMDistribution(m);
  init_Absorption(m);
  init_STTSLReaction(m);
  init_STIncoherentInelastic(m);
  init_STCoherentElastic(m);
  init_STInoherentElasticACE(m);
  init_STInoherentElastic(m);
  init_STThermalScatteringLaw(m);
  init_URRPtable(m);
  init_PCTable(m);
  init_EnergyAngleTable(m);
  init_TabularEnergyAngle(m);
  init_EnergyGrid(m);
  init_CrossSection(m);
  init_ReactionBase(m);
  init_STReaction(m);
  init_DelayedFamily(m);
  init_Fission(m);
  init_STNeutron(m);
  init_PRNG(m);
  init_ZAID(m);
  init_Element(m);
  init_Isotope(m);
  init_Nuclide(m);
  init_NDLibrary(m);
  init_MCNPLibrary(m);
  init_SerpentLibrary(m);

  m.attr("__author__") = "Hunter Belanger";
  m.attr("__copyright__") = "Copyright 2021-2023, Hunter Belanger";
  m.attr("__license__") = "GPL-3.0-or-later";
  m.attr("__maintainer__") = "Hunter Belanger";
  m.attr("__email__") = "hunter.belanger@gmail.com";
  m.attr("__version__") = pndl::VERSION_STRING;
}
