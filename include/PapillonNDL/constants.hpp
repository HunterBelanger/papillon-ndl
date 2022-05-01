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
#ifndef PAPILLON_NDL_CONSTANTS_H
#define PAPILLON_NDL_CONSTANTS_H

namespace pndl {

constexpr double SEC_TO_SHAKE = 1.E-8;
constexpr double SHAKE_TO_SEC = 1.E8;
constexpr double EV_TO_K = 1.160451812E4;
constexpr double MEV_TO_EV = 1.E6;
constexpr double EV_TO_MEV = 1.E-6;
constexpr double KB = 8.61733034E-11;  // [MeV / K]
constexpr double PI = 3.14159265358979323846264338327950288;
constexpr double SQRT_PI = 1.77245385090551602729816748334114518;

}

#endif  // PAPILLON_NDL_CONSTANTS_H
