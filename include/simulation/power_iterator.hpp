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
#ifndef MG_POWER_ITERATOR_H
#define MG_POWER_ITERATOR_H

#include <simulation/simulation.hpp>

class PowerIterator : public Simulation {
 public:
  PowerIterator(std::shared_ptr<Tallies> i_t, std::shared_ptr<Settings> i_s,
                std::shared_ptr<Transporter> i_tr,
                std::vector<std::shared_ptr<Source>> src)
      : Simulation{i_t, i_s, i_tr, src}, bank{} {};
  ~PowerIterator() = default;

  void initialize();

  void run();

  void premature_kill();

 private:
  std::vector<Particle> bank;

  void generation_output(int gen, int Nn, int Np, int Npos, int Nneg, int Wn,
                         int Wp, int Wpos, int Wneg);

  void write_source(std::vector<Particle>& bank, std::string source_fname);

  std::vector<std::string> split(std::string line);

  std::string remove(char c, std::string line);

  void perform_regional_cancelation(std::vector<Particle>& next_gen);

  void print_header();

};  // PowerIterator

#endif  // MG_POWER_ITERATOR_H
