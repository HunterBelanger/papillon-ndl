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
#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <fstream>

#include "constants.hpp"

namespace pndl {
// Forward declaration of split_line function
static std::vector<std::string> split_line(std::string line);

ACE::ACE(std::string fname)
    : zaid_(),
      temperature_(),
      awr_(),
      fissile_(),
      izaw_(),
      nxs_(),
      jxs_(),
      xss_(),
      esz_(),
      nu_(),
      mtr_(),
      lqr_(),
      tyr_(),
      lsig_(),
      sig_(),
      lan_(),
      an_(),
      ldlw_(),
      dlw_(),
      dnedl_(),
      dned_(),
      dnu_(),
      bdd_() {
  // Open ACE file
  std::ifstream file(fname);

  std::string line;
  std::getline(file, line);

  // Check first line to determine header type
  bool legacy_header = true;
  if (line[0] == '2' && line[1] == '.') legacy_header = false;

  std::vector<std::string> split;

  // Parse header
  if (legacy_header) {
    split = split_line(line);
    awr_ = std::stod(split[1]);
    temperature_ = std::stod(split[2]) * MEV_TO_EV * EV_TO_K;

    // Skip next line
    std::getline(file, line);
  } else {
    // Read next line
    std::getline(file, line);
    split = split_line(line);
    awr_ = std::stod(split[0]);
    temperature_ = std::stod(split[1]) * MEV_TO_EV * EV_TO_K;
    int n_skip = std::stoi(split[3]);

    // Skip comment lines
    for (int i = 0; i < n_skip; i++) std::getline(file, line);
  }

  // Parse IZAW
  for (int i = 0; i < 16; i++) {
    int32_t i_zaid;
    double i_awr;
    file >> i_zaid;
    file >> i_awr;
    izaw_[i] = {i_zaid, i_awr};
  }

  // Parse NXS
  for (int i = 0; i < 16; i++) {
    file >> nxs_[i];
  }

  // Parse JXS
  for (int i = 0; i < 32; i++) {
    file >> jxs_[i];
  }

  // Parse XSS
  xss_.resize(nxs_[0]);
  int i = 0;
  while(!file.eof()) {
    file >> xss_[i];
    i++;
  }

  if(i-1 != nxs_[0]) {
    std::string mssg = "ACE::ACE: Found incorrect number of entries in XSS array while reading\n";
    mssg +=            "the \"" + fname + "\" ACE file.\n";
    mssg +=            "This is likely due to a numerical entry which is missing the \"E\".\n";
    mssg +=            "Please correct the ACE file.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  zaid_ = static_cast<uint32_t>(nxs_[1]);

  if (jxs_[1] > 0) fissile_ = true;

  // Set locator constants
  esz_ = jxs_[0] - 1;
  nu_ = jxs_[1] - 1;
  mtr_ = jxs_[2] - 1;
  lqr_ = jxs_[3] - 1;
  tyr_ = jxs_[4] - 1;
  lsig_ = jxs_[5] - 1;
  sig_ = jxs_[6] - 1;
  lan_ = jxs_[7] - 1;
  an_ = jxs_[8] - 1;
  ldlw_ = jxs_[9] - 1;
  dlw_ = jxs_[10] - 1;
  dnu_ = jxs_[23] - 1;
  bdd_ = jxs_[24] - 1;
  dnedl_ = jxs_[25] - 1;
  dned_ = jxs_[26] - 1;
}

int32_t ACE::zaid() const { return zaid_; }
double ACE::temperature() const { return temperature_; }
double ACE::awr() const { return awr_; }
bool ACE::fissile() const { return fissile_; }

std::pair<int32_t, double> ACE::izaw(size_t i) const { return izaw_[i]; }
int32_t ACE::nxs(size_t i) const { return nxs_[i]; }
int32_t ACE::jxs(size_t i) const { return jxs_[i]; }
double ACE::xss(size_t i) const { return xss_[i]; }

std::vector<std::pair<int32_t, double>> ACE::izaw(size_t i, size_t len) const {
  return {izaw_.begin() + i, izaw_.begin() + i + len};
}

std::vector<int32_t> ACE::nxs(size_t i, size_t len) const {
  return {nxs_.begin() + i, nxs_.begin() + i + len};
}

std::vector<int32_t> ACE::jxs(size_t i, size_t len) const {
  return {jxs_.begin() + i, jxs_.begin() + i + len};
}

std::vector<double> ACE::xss(size_t i, size_t len) const {
  return {xss_.begin() + i, xss_.begin() + i + len};
}

const double* ACE::xss_data() const { return xss_.data(); }

int32_t ACE::ESZ() const { return esz_; }
int32_t ACE::NU() const { return nu_; }
int32_t ACE::MTR() const { return mtr_; }
int32_t ACE::LQR() const { return lqr_; }
int32_t ACE::TYR() const { return tyr_; }
int32_t ACE::LSIG() const { return lsig_; }
int32_t ACE::SIG() const { return sig_; }
int32_t ACE::LAND() const { return lan_; }
int32_t ACE::AND() const { return an_; }
int32_t ACE::LDLW() const { return ldlw_; }
int32_t ACE::DLW() const { return dlw_; }
int32_t ACE::DNEDL() const { return dnedl_; }
int32_t ACE::DNED() const { return dned_; }
int32_t ACE::DNU() const { return dnu_; }
int32_t ACE::BDD() const { return bdd_; }

static std::vector<std::string> split_line(std::string line) {
  std::vector<std::string> out;

  std::string tmp = "";
  for (size_t i = 0; i < line.size(); i++) {
    if (line[i] != ' ')
      tmp += line[i];
    else {
      if (tmp.size() > 0) {
        out.push_back(tmp);
        tmp = "";
      }
    }
  }

  if (tmp.size() > 0) out.push_back(tmp);

  return out;
}

}  // namespace pndl
