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
