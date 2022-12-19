/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2022, Hunter Belanger
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

/**
 * @file
 * @author Hunter Belanger
 */

#include "coherent_elastic.hpp"
#include "incoherent_elastic.hpp"
#include "incoherent_inelastic.hpp"
#include "linearize.hpp"
#include "ace.hpp"


#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <variant>


#include <ENDFtk.hpp>
using namespace njoy::ENDFtk;
#include <Log.hpp>
using namespace njoy;
#include <boost/hana.hpp>  // Needed for the _c literal for constructing mt4
#include <docopt.h>

#define VERSION_STRING "0.1.0"

static const std::string version =
  "Panglos : A Thermal Scattering Law Processor\n"
  "Version " VERSION_STRING "\n\n"

  "Copyright (C) 2022 Hunter Belanger.\n"
  "Released under the terms and conditions of the GPLv3 license.\n"
  "Written by Hunter Belanger.\n";

static const std::string usage =
  "Usage:\n"
  "  panglos process [--pedantic] <fname> <mat> <temp> <zaid> <comments> <acefname>\n"
  "  panglos temps <fname> <mat>\n"
  "  panglos (-h | --help)\n"
  "  panglos (-v | --version)\n\n"
  
  "Options:\n"
  "  -p --pedantic  Perform pedantic checks on distribution linearization\n"
  "  -h --help      Show this help message\n"
  "  -v --version   Show version number\n";

static const std::string help = version + '\n' + usage;

void read_elastic(const file::Type<7>& mf7,
                  std::unique_ptr<CoherentElastic>& ce,
                  std::unique_ptr<IncoherentElastic>& ie) {
  if (mf7.hasSection(2)) {
    section::Type<7, 2> mt2 = mf7.section(2_c);
    const auto& scatter_law = mt2.scatteringLaw();

    if (std::holds_alternative<section::Type<7, 2>::CoherentElastic>(
            scatter_law)) {
      ce = std::make_unique<CoherentElastic>(
          std::get<section::Type<7, 2>::CoherentElastic>(scatter_law));
    } else if (std::holds_alternative<section::Type<7, 2>::IncoherentElastic>(
                   scatter_law)) {
      ie = std::make_unique<IncoherentElastic>(
          std::get<section::Type<7, 2>::IncoherentElastic>(scatter_law));
    } else {
      const section::Type<7, 2>::MixedElastic& me =
          std::get<section::Type<7, 2>::MixedElastic>(scatter_law);
      ce = std::make_unique<CoherentElastic>(me.coherent());
      ie = std::make_unique<IncoherentElastic>(me.incoherent());
    }
  }
}

int main(const int argc, const char** argv) {
  // Initialize docopt
  std::map<std::string, docopt::value> args =
      docopt::docopt(help, {argv + 1, argv + argc}, true, version);

  // Read parameters
  const std::string fname = args["<fname>"].asString();
  const int MAT = static_cast<int>(args["<mat>"].asLong());
  const bool get_temps = args["temps"].asBool();
  const bool pedantic = args["--pedantic"].asBool();
  

  // Write run options
  Log::info("");
  Log::info("Panglos : A Thermal Scattering Law Processor");
  Log::info("-----------------------------------------------------------");
  Log::info("Copyright (C) 2022 Hunter Belanger");
  Log::info("Released under the terms and conditions of the GPLv3.");
  Log::info("");
  Log::info("File Name:   {}", fname);
  Log::info("MAT:         {}", MAT); 

  //=============================================================================
  // Check if we were asked to just list the temperatures
  if (get_temps) {
    // Read ENDF file, and get MF7
    tree::Tape<std::string> endf = tree::fromFile(fname);
    file::Type<7> mf7 = endf.material(MAT).front().file(7).parse<7>();

    // First read incoherent inelastic, which must be present
    section::Type<7, 4> mt4 = mf7.section(4_c);
    IncoherentInelastic ii(mt4);

    std::stringstream temps;
    temps << "[";
    const std::size_t NT = ii.temperatures().size();
    for (std::size_t i = 0; i < NT; i++) {
      temps << std::fixed << std::setprecision(1);
      temps << ii.temperatures()[i];
      if (i != NT-1) {
        temps << ", ";
      }
    }
    temps << "]";

    Log::info("Provided Temperatures: {}", temps.str());
    return 0;
  }

  //=============================================================================
  const double T = std::stod(args["<temp>"].asString());
  std::string zaid = args["<zaid>"].asString();
  std::string comments = args["<comments>"].asString();
  std::string acefname = args["<acefname>"].asString();

  Log::info("Temperature: {}", T);
  Log::info("ZAID:        {}", zaid);
  Log::info("Comments:    {}", comments);
  Log::info("ACE File:    {}", acefname);
  if (pedantic) {
    Log::info("Pedantic:    True");
  } else {
    Log::info("Pedantic:    False");
  }
  Log::info("");

  if (zaid.size() < 10) {
    while (zaid.size() < 10) {
      zaid.insert(zaid.begin(), ' ');
    }
  } else if (zaid.size() > 10) {
    zaid.resize(10, ' ');
  }
  comments.resize(70, ' ');

  // Read ENDF file, and get MF7
  tree::Tape<std::string> endf = tree::fromFile(fname);
  file::Type<7> mf7 = endf.material(MAT).front().file(7).parse<7>();

  // First read incoherent inelastic, which must be present
  section::Type<7, 4> mt4 = mf7.section(4_c);
  IncoherentInelastic ii(mt4);

  // Check for our temperature
  std::size_t Ti = ii.temperatures().size();
  for (std::size_t i = 0; i < ii.temperatures().size(); i++) {
    const double abs_temp_diff = std::abs(T - ii.temperatures()[i]);
    if (abs_temp_diff < 1.) {
      Ti = i;
      break;
    }
  }
  if (Ti == ii.temperatures().size()) {
    Log::error("Could not find a tabulated scattering law for {} K.", T);
    return 1;
  }
  Log::info("Processing at temperature {} K.", ii.temperatures()[Ti]);
  Log::info("");

  // If elastic components are present, read those too
  std::unique_ptr<CoherentElastic> ce = nullptr;
  std::unique_ptr<IncoherentElastic> ie = nullptr;
  read_elastic(mf7, ce, ie);
  
  // Linearize the Incoherent Inelastic xs and distribution
  LinearizedIncoherentInelastic lii = linearize_ii(ii, Ti, pedantic);

  // Write data to ACE file
  Log::info("Writing ACE file.");
  write_to_ace(lii, ie, ce, zaid, ii.awr(), T, comments, MAT, acefname);
  Log::info("");

  Log::info("TSL Processing Complete !");
  Log::info("");

  return 0;
}
