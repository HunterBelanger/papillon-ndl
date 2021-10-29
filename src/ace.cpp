/*
 * Papillon Nuclear Data Library
 * Copyright 2021, Hunter Belanger
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
#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/constants.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <filesystem>
#include <fstream>
#include <ios>
#include <string>


namespace pndl {

// Forward declaration of split_line function
static std::vector<std::string> split_line(std::string line);

ACE::ACE(std::string fname, Type type)
    : zaid_(),
      temperature_(),
      awr_(),
      fissile_(),
      fname_(fname),
      zaid_txt(10, ' '),
      date_(10, ' '),
      comment_(70, ' '),
      mat_(10, ' '),
      izaw_(),
      nxs_(),
      jxs_(),
      xss_() {
  // Make sure file exists
  if (!std::filesystem::exists(fname)) {
    std::string mssg = "File \"" + fname + "\" does not exist.";
    throw PNDLException(mssg);
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
  // Check first line to determine header type
  bool legacy_header = true;

  char c1 = static_cast<char>(file.get());
  char c2 = static_cast<char>(file.peek());
  if (c1 == '2' && c2 == '.') legacy_header = false;
  file.unget();

  // Parse header
  std::string awr_txt(12, ' ');
  std::string temp_txt(12, ' ');
  if (legacy_header) {
    file.read(zaid_txt.data(), 10);
    file.read(awr_txt.data(), 12);
    file.read(temp_txt.data(), 12);
    awr_ = std::stod(awr_txt);
    temperature_ = std::stod(temp_txt) * MEV_TO_EV * EV_TO_K;

    // Skip blank char
    file.ignore(1);

    // Read date
    file.read(date_.data(), 10);

    // Ignore the newline chars
    if (file.peek() == '\n' || file.peek() == '\r') file.ignore(1);
    if (file.peek() == '\n' || file.peek() == '\r') file.ignore(1);

    // Read comment
    file.read(comment_.data(), 70);

    // Read mat id
    file.read(mat_.data(), 10);
  } else {
    std::string line;
    std::getline(file, line);
    // Read next line
    std::getline(file, line);
    std::vector<std::string> split = split_line(line);
    awr_ = std::stod(split[0]);
    temperature_ = std::stod(split[1]) * MEV_TO_EV * EV_TO_K;
    int n_skip = std::stoi(split[3]);

    if (n_skip == 2) {
      // These are the legacy header. Read them
      // Read zaid text
      file.read(zaid_txt.data(), 10);

      // Skip duplicate awr and temp and space
      file.ignore(25);

      // Read date
      file.read(date_.data(), 10);

      // Ignore the newline chars
      if (file.peek() == '\n' || file.peek() == '\r') file.ignore(1);
      if (file.peek() == '\n' || file.peek() == '\r') file.ignore(1);

      // Read comment
      file.read(comment_.data(), 70);

      // Read mat id
      file.read(mat_.data(), 10);
    } else {
      // Skip comment lines
      for (int i = 0; i < n_skip; i++) std::getline(file, line);
    }
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
  while (!file.eof() && i < nxs_[0]) {
    file >> xss_[i];
    i++;
  }

  if (i != nxs_[0]) {
    std::string mssg =
        "Found incorrect number of entries in XSS array while reading the \"" +
        fname_ +
        "\" ACE file. This is likely due to a numerical entry which is missing "
        "the \"E\". Please correct the ACE file.";
    throw PNDLException(mssg);
  }

  zaid_ = static_cast<uint32_t>(nxs_[1]);

  if (jxs_[1] > 0) fissile_ = true;
}

void ACE::read_binary(std::ifstream& file) {
  // Skip first record length
  file.ignore(4);

  // Skip zaid
  file.read(zaid_txt.data(), 10);

  // Read the AWR
  file.read(reinterpret_cast<char*>(&awr_), sizeof(double));

  // Read the temperatuer
  file.read(reinterpret_cast<char*>(&temperature_), sizeof(double));
  temperature_ *= MEV_TO_EV * EV_TO_K;

  // Read date
  file.read(date_.data(), 10);

  // Read comment
  file.read(comment_.data(), 70);

  // Read mat
  file.read(mat_.data(), 10);

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

void ACE::save_binary(std::string& fname) {
  std::ofstream file(fname, std::ios_base::binary);

  // Write first record length which is size of all
  // of the ACE header
  uint32_t rlen = 100 + 64 * sizeof(int32_t) + 18 * sizeof(double);
  file.write(reinterpret_cast<char*>(&rlen), 4);

  // Write zaid
  file.write(zaid_txt.data(), 10);

  // Write the AWR
  file.write(reinterpret_cast<char*>(&awr_), sizeof(double));

  // Write the temperatuer
  temperature_ /= MEV_TO_EV * EV_TO_K;
  file.write(reinterpret_cast<char*>(&temperature_), sizeof(double));

  // Write date, comment, and mat
  file.write(date_.data(), 10);
  file.write(comment_.data(), 70);
  file.write(mat_.data(), 10);

  // Write IZAW
  for (std::size_t i = 0; i < 16; i++) {
    file.write(reinterpret_cast<char*>(&izaw_[i].first), sizeof(int32_t));
    file.write(reinterpret_cast<char*>(&izaw_[i].second), sizeof(double));
  }

  // Write NXS
  for (std::size_t i = 0; i < 16; i++) {
    file.write(reinterpret_cast<char*>(&nxs_[i]), sizeof(int32_t));
  }

  // Write JXS
  for (std::size_t i = 0; i < 32; i++) {
    file.write(reinterpret_cast<char*>(&jxs_[i]), sizeof(int32_t));
  }

  // Write end record length
  file.write(reinterpret_cast<char*>(&rlen), 4);

  // Write XSS
  const uint32_t ner = 512;
  std::size_t ll = 0;
  std::size_t nn = xss_.size();
  while (nn > 0) {
    std::size_t n = nn;

    if (n > ner) n = ner;

    rlen = static_cast<uint32_t>(n * sizeof(double));

    // Write first len header
    file.write(reinterpret_cast<char*>(&rlen), 4);

    for (std::size_t j = ll; j < ll + n; j++) {
      file.write(reinterpret_cast<char*>(&xss_[j]), sizeof(double));
    }
    ll += n;
    nn -= n;

    // Write second len header
    file.write(reinterpret_cast<char*>(&rlen), 4);
  }

  file.close();
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
