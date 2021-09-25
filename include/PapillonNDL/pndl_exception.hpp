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
#ifndef PAPILLON_NDL_EXCEPTION_H
#define PAPILLON_NDL_EXCEPTION_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <exception>
#include <source_location>
#include <string>

namespace pndl {

/**
 * @brief Class used for exceptions by the library.
 */
class PNDLException : public std::exception {
 public:
  PNDLException() : message("\n") {}
  /**
   * @param mssg Error message.
   * @param file Location where the error occurred;
   */
  PNDLException(const std::string& mssg,
                std::source_location location = std::source_location::current())
      : message("\n") {
    add_to_error_message(mssg, location);
  }
  ~PNDLException() = default;

  /**
   * @brief Adds details to the exception message as it is passed up the stack.
   * @param mssg Error message.
   * @param file Location where the error was thrown.
   */
  void add_to_exception(
      const std::string& mssg,
      std::source_location location = std::source_location::current()) {
    add_to_error_message(mssg, location);
  }

  const char* what() const noexcept override { return message.c_str(); }

 private:
  std::string message;

  void add_to_error_message(const std::string& mssg,
                            const std::source_location& location) {
    // Go through original string and determine line breaks
    std::string mssg_tmp = mssg;
    int nbreaks = static_cast<int>(mssg_tmp.size()) / 80;
    if ((mssg_tmp.size() % 80) == 0) nbreaks--;
    if (nbreaks > 0) {
      for (size_t b = 0; b < static_cast<size_t>(nbreaks); b++) {
        // Get index of break.
        size_t ind = 80 * (b + 1);

        // Work backwards to first space
        while (mssg_tmp[ind] != ' ') {
          if (ind > 0)
            ind--;
          else {
            ind = 80 * (b + 1);
            break;
          }
        }

        mssg_tmp[ind] = '\n';
      }
    }

    std::string tmp = "\n";
    tmp +=
        " #--------------------------------------------------------------------"
        "-------------\n";
    tmp += " # File: " + std::string(location.file_name()) + "\n";
    tmp += " # Function: " + std::string(location.function_name()) + "\n";
    tmp += " # Line: " + std::to_string(location.line()) + "\n";
    tmp += " # \n";
    tmp += " # ";

    for (const auto& c : mssg_tmp) {
      if (c == '\n') {
        tmp += "\n # ";
      } else {
        tmp += c;
      }
    }
    tmp += "\n";
    tmp +=
        " #--------------------------------------------------------------------"
        "-------------";

    message = tmp + message;
  }
};
}  // namespace pndl

#endif
