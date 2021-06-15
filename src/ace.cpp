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
#include <filesystem>
#include <fstream>
#include <ios>

#include "constants.hpp"

namespace pndl {
// Forward declaration of split_line function
static std::vector<std::string> split_line(std::string line);

ACE::ACE(std::string fname, Type type)
    : zaid_(),
      temperature_(),
      awr_(),
      fissile_(),
      fname_(fname),
      izaw_(),
      nxs_(),
      jxs_(),
      xss_() {
  // Make sure file exists
  if (!std::filesystem::exists(fname)) {
    std::string mssg = "ACE::ACE: File \"" + fname + "\" does not exist.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  // Open ACE file
  std::ios_base::openmode mode = std::ios_base::in;
  if (type == Type::BINARY) mode = std::ios_base::binary;
  std::ifstream file(fname, mode);

  switch (type) {
    case Type::ASCII:
      read_ascii(file);
      break;

    case Type::BINARY:
      read_binary(file);
      break;
  }

  file.close();
}

void ACE::read_ascii(std::ifstream& file) {
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
  while (!file.eof() && i <= nxs_[0]) {
    file >> xss_[i];
    i++;
  }

  if (i - 1 != nxs_[0]) {
    std::string mssg =
        "ACE::read_ascii: Found incorrect number of entries in XSS array while "
        "reading the \"" +
        fname_ +
        "\" ACE file. This is likely due to a numerical entry which is missing "
        "the \"E\". Please correct the ACE file.";
    throw PNDLException(mssg, __FILE__, __LINE__);
  }

  zaid_ = static_cast<uint32_t>(nxs_[1]);

  if (jxs_[1] > 0) fissile_ = true;
}

void ACE::read_binary(std::ifstream& file) {
  // Skip first record length
  file.ignore(4);

  // Skip zaid
  file.ignore(10);

  // Read the AWR
  file.read(reinterpret_cast<char*>(&awr_), sizeof(double));

  // Read the temperatuer
  file.read(reinterpret_cast<char*>(&temperature_), sizeof(double));
  temperature_ *= MEV_TO_EV * EV_TO_K;

  // Skip date, comment, and mat
  file.ignore(90);

  // Parse IZAW
  for (int i = 0; i < 16; i++) {
    int32_t i_zaid;
    double i_awr;

    file.read(reinterpret_cast<char*>(&i_zaid), sizeof(int32_t));
    file.read(reinterpret_cast<char*>(&i_awr), sizeof(double));
    izaw_[i] = {i_zaid, i_awr};
  }

  // Parse NXS
  for (int i = 0; i < 16; i++) {
    file.read(reinterpret_cast<char*>(&nxs_[i]), sizeof(int32_t));
  }

  // Parse JXS
  for (int i = 0; i < 32; i++) {
    file.read(reinterpret_cast<char*>(&jxs_[i]), sizeof(int32_t));
  }

  // Skip end record length
  file.ignore(4);

  // Parse XSS
  xss_.resize(nxs_[0]);
  uint32_t rlen;
  std::size_t i = 0;
  while (!file.eof() && i < xss_.size()) {
    // Get the record entry length
    file.read(reinterpret_cast<char*>(&rlen), 4);

    file.read(reinterpret_cast<char*>(&xss_[i]), rlen);
    i += rlen / sizeof(double);

    // Skip last len entry
    file.ignore(4);
  }

  zaid_ = static_cast<uint32_t>(nxs_[1]);

  if (jxs_[1] > 0) fissile_ = true;
}

std::vector<std::pair<int32_t, double>> ACE::izaw(std::size_t i,
                                                  std::size_t len) const {
  return {izaw_.begin() + i, izaw_.begin() + i + len};
}

std::vector<int32_t> ACE::nxs(std::size_t i, std::size_t len) const {
  return {nxs_.begin() + i, nxs_.begin() + i + len};
}

std::vector<int32_t> ACE::jxs(std::size_t i, std::size_t len) const {
  return {jxs_.begin() + i, jxs_.begin() + i + len};
}

std::vector<double> ACE::xss(std::size_t i, std::size_t len) const {
  return {xss_.begin() + i, xss_.begin() + i + len};
}

const double* ACE::xss_data() const { return xss_.data(); }

static std::vector<std::string> split_line(std::string line) {
  std::vector<std::string> out;

  std::string tmp = "";
  for (std::size_t i = 0; i < line.size(); i++) {
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
