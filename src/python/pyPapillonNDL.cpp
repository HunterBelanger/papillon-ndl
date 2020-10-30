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
#include <pybind11/pybind11.h>

namespace py = pybind11;

extern void init_ACE(py::module&);
extern void init_Interpolation(py::module&);
extern void init_Frame(py::module&);
extern void init_Function1D(py::module&);
extern void init_Constant(py::module&);
extern void init_Polynomial1D(py::module&);
extern void init_Tabulated1D(py::module&);
extern void init_Region1D(py::module&);
extern void init_MultiRegion1D(py::module&);
extern void init_AngleLaw(py::module&);
extern void init_Isotropic(py::module&);
extern void init_EquiprobableAngleBins(py::module&);
extern void init_AngleTable(py::module&);
extern void init_AngleDistribution(py::module&);
extern void init_EnergyLaw(py::module&);
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
extern void init_NBody(py::module&);
extern void init_KalbachTable(py::module&);
extern void init_Kalbach(py::module&);
extern void init_PCTable(py::module&);
extern void init_EnergyAngleTable(py::module&);
extern void init_TabularEnergyAngle(py::module&);
extern void init_SharedSpanFloat(py::module&);
extern void init_EnergyGrid(py::module&);
extern void init_CrossSection(py::module&);
extern void init_Reaction(py::module&);
extern void init_Nuclide(py::module&);

PYBIND11_MODULE(pyPapillonNDL, m){
  init_ACE(m); 
  init_Interpolation(m);
  init_Frame(m);
  init_Function1D(m);
  init_Constant(m);
  init_Polynomial1D(m);
  init_Tabulated1D(m);
  init_Region1D(m);
  init_MultiRegion1D(m);
  init_AngleLaw(m);
  init_Isotropic(m);
  init_EquiprobableAngleBins(m);
  init_AngleTable(m);
  init_AngleDistribution(m);
  init_EnergyLaw(m);
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
  init_NBody(m);
  init_KalbachTable(m);
  init_Kalbach(m);
  init_PCTable(m);
  init_EnergyAngleTable(m);
  init_TabularEnergyAngle(m);
  init_SharedSpanFloat(m);
  init_EnergyGrid(m);
  init_CrossSection(m);
  init_Reaction(m);
  init_Nuclide(m);
}
