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
#include <simulation/tallies.hpp>
#include <vector>

Tallies::Tallies(double tot_wgt)
    :

      total_weight(tot_wgt),
      k_col_score(0.),
      k_tot_score(0.),
      leak_score(0.),
      k_eff(1.0),
      k_tot(1.0),
      k_avg(0.0),
      k_var(0.0),
      leak(0.0),
      leak_avg(0.0),
      leak_var(0.0) {}

void Tallies::set_flux_tally(std::unique_ptr<FluxTally> ftally) {
  flux_tally = std::move(ftally);
}

void Tallies::set_power_tally(std::unique_ptr<PowerTally> ptally) {
  power_tally = std::move(ptally);
}

void Tallies::score_collision(const Particle& p,
                              const std::shared_ptr<Material>& mat,
                              bool converged) {
  double Et = mat->Et(p.r(), p.E());
  double Ef = mat->Ef(p.r(), p.E());
  double nu = mat->nu(p.r(), p.E());

  double scr = p.wgt() * nu * Ef / Et;

#pragma omp atomic
  k_col_score += scr;
#pragma omp atomic
  k_tot_score += std::abs(scr);

  // Only do spacial tallies if converged
  if (converged) {
    if (flux_tally) flux_tally->score_flux(p, mat);

    if (power_tally) power_tally->score_power(p, mat);
  }
}

void Tallies::score_leak(double scr) {
#pragma omp atomic
  leak_score += scr;
}

void Tallies::clear_generation(bool converged) {
  k_col_score = 0.;
  k_tot_score = 0.;
  leak_score = 0.;

  if (converged) {
    if (flux_tally) flux_tally->clear_generation();

    if (power_tally) power_tally->clear_generation();
  }
}

void Tallies::calc_gen_values() {
  k_eff = k_col_score / total_weight;
  k_tot = k_tot_score / total_weight;
  leak = leak_score / total_weight;
}

void Tallies::record_generation() {
  gen++;

  double old_k_avg = k_avg;
  k_avg = k_avg + (k_eff - k_avg) / static_cast<double>(gen);
  k_var = k_var + ((k_eff - old_k_avg) * (k_eff - k_avg) - (k_var)) /
                      static_cast<double>(gen);

  double old_leak_avg = leak_avg;
  leak_avg = leak_avg + (leak - leak_avg) / static_cast<double>(gen);
  leak_var =
      leak_var + ((leak - old_leak_avg) * (leak - leak_avg) - (leak_var)) /
                     static_cast<double>(gen);

  if (flux_tally) flux_tally->record_generation(static_cast<double>(gen));

  if (power_tally) power_tally->record_generation(static_cast<double>(gen));
}

void Tallies::write_flux(std::string flux_fname) {
  if (flux_tally) flux_tally->write_flux(flux_fname);
}

void Tallies::write_power(std::string power_fname) {
  if (power_tally) power_tally->write_power(power_fname);
}
