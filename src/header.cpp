/*=============================================================================*
 * Copyright (C) 2021, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
 *
 * Contributeur : Hunter Belanger (hunter.belanger@cea.fr)
 *
 * Ce logiciel est un programme informatique servant à faire des comparaisons
 * entre les méthodes de transport qui sont capable de traiter les milieux
 * continus avec la méthode Monte Carlo. Il résoud l'équation de Boltzmann
 * pour les particules neutres, à une vitesse et dans une dimension.
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
 *============================================================================*/
#include <utils/constants.hpp>
#include <utils/header.hpp>

const std::string logo =
    "\n"
    "                         ███╗   ███╗ ██████╗ ███╗   ███╗ ██████╗\n"
    "                         ████╗ ████║██╔════╝ ████╗ ████║██╔════╝\n"
    "                         ██╔████╔██║██║  ███╗██╔████╔██║██║     \n"
    "                         ██║╚██╔╝██║██║   ██║██║╚██╔╝██║██║     \n"
    "                         ██║ ╚═╝ ██║╚██████╔╝██║ ╚═╝ ██║╚██████╗\n"
    "                         ╚═╝     ╚═╝ ╚═════╝ ╚═╝     ╚═╝ ╚═════╝\n\n";

const std::string header =
    "                         Multi-Group Monte Carlo Transport Code\n\n";

const std::string info =
    " Copyright (C) 2021, Commissariat à l'énergie atomique et aux énergies "
    "alternatives (CEA)\n"
    " Released under the terms and conditions of the CeCILL-v2.1 license\n"
    " Written by Hunter Belanger\n";

const std::string help =
    " Usage:\n"
#ifdef _OPENMP
    "   mgmc (--input FILE) [--threads NUM --output FILE]\n"
#else
    "   mgmc (--input FILE) [--output FILE]\n"
#endif
    "   mgmc (--input FILE --plot) [--threads NUM]\n"
    "   mgmc (-h | --help)\n"
    "   mgmc (-v | --version)\n\n"

    " Options:\n"
    "   -h --help         Show this help message\n"
    "   -v --version      Show version number\n"
    "   -i --input FILE   Set input file\n"
#ifdef _OPENMP
    "   -t --threads NUM  Set number of OpenMP threads\n"
#endif
    "   -o --output FILE  Set output file\n"
    "   -p --plot         Generate plots from input file\n";

const std::string version_string = " MGMC " MGMC_VERSION_STRING
                                   " "
#if defined(DEV_VERSION)
                                   "(Development)"
#endif
                                   "\n"
                                   " Git Hash " MGMC_GIT_HASH
                                   "\n"
#if defined(_OPENMP)
                                   " Compiled with OpenMP.\n"
#endif
#if defined(MGMC_USE_MPI)
                                   " Compiled with MPI.\n"
#endif
                                   "";
