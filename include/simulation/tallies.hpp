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
#ifndef MG_TALLIES_H
#define MG_TALLIES_H

#include <materials/material.hpp>
#include <simulation/flux_tally.hpp>
#include <simulation/particle.hpp>
#include <simulation/power_tally.hpp>
#include <string>

class Tallies {
 public:
  Tallies(double tot_wgt);
  ~Tallies() = default;

  void set_flux_tally(std::unique_ptr<FluxTally> ftally);

  void set_power_tally(std::unique_ptr<PowerTally> ptally);

  void score_collision(const Particle& p, const std::shared_ptr<Material>& mat,
                       bool converged);

  void score_leak(double scr);

  void clear_generation(bool converged);

  void calc_gen_values();

  void record_generation();

  double keff() const { return k_eff; }
  double ktot() const { return k_tot; }
  double kavg() const { return k_avg; }
  double kerr() const { return std::sqrt(k_var / (static_cast<double>(gen))); }

  void write_flux(std::string flux_fname);

  void write_power(std::string power_fname);

  void set_total_weight(double tot_wgt) { total_weight = tot_wgt; }

 private:
  int gen = 0;

  double total_weight;

  double k_col_score;
  double k_tot_score;
  double leak_score;

  double k_eff, k_tot;
  double k_avg, k_var;
  double leak, leak_avg, leak_var;

  std::unique_ptr<FluxTally> flux_tally = nullptr;
  std::unique_ptr<PowerTally> power_tally = nullptr;

};  // Tallies

#endif  // MG_TALLIES_H
