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
#ifndef PANGLOS_CONSTANTS_H
#define PANGLOS_CONSTANTS_H

/**
 * @file
 * @author Hunter Belanger
 */

constexpr double KB = 8.617333262E-5;  // Boltzmann Constant [eV / K]
constexpr double EV_TO_MEV = 1.0E-6; // 10^(-6) MeV / eV = 1
constexpr double K_TO_MEV = KB * EV_TO_MEV;
constexpr double TROOM = 0.0253 / KB;  // Room Temperature [K]
constexpr double PI = 3.1415926535897932384626433832795028841971694;
constexpr double SCT_CUTOFF = 1.921947727823849e-98;  // exp(-225)

#endif
