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
#ifndef TRANSPORTER_H
#define TRANSPORTER_H

#include <geometry/geometry.hpp>
#include <materials/nuclide.hpp>
#include <optional>
#include <simulation/noise_source.hpp>
#include <simulation/particle.hpp>
#include <simulation/tallies.hpp>
#include <utils/constants.hpp>
#include <utils/rng.hpp>
#include <utils/settings.hpp>

class Transporter {
 public:
  Transporter(std::shared_ptr<Tallies> i_t) : tallies{i_t} {};
  virtual ~Transporter() = default;

  virtual std::vector<BankedParticle> transport(
      std::vector<Particle> &bank, bool noise = false,
      std::vector<BankedParticle> *noise_bank = nullptr,
      std::vector<std::shared_ptr<NoiseSource>> *noise_sources = nullptr) = 0;

 protected:
  std::shared_ptr<Tallies> tallies;

  struct ThreadLocalScores {
    double k_col_score = 0.;
    double k_abs_score = 0.;
    double k_trk_score = 0.;
    double k_tot_score = 0.;
    double leakage_score = 0.;
  };

  void russian_roulette(Particle &p);

  void collision(
      Particle &p, MaterialHelper &mat, ThreadLocalScores &thread_scores,
      bool noise = false,
      std::vector<std::shared_ptr<NoiseSource>> *noise_sources = nullptr);

  // Function to get the speed of a particle in cm/s
  double speed(double E, std::size_t i) {
    if (settings::energy_mode == settings::EnergyMode::MG) {
      return settings::mg_speeds[i];
    } else {
      double m = N_MASS_EV / (C_CM_S * C_CM_S);  // Mass in [eV * s^2 / cm^2]
      return std::sqrt(2. * E / m);              // Speed in [cm / s]
    }
  }

  void make_noise_copy(Particle &p, const MicroXSs &microxs) const;

  void make_fission_neutrons(Particle &p, const MicroXSs &microxs,
                             const Nuclide &nuclide, bool noise) const;

  void sample_branchless_noise(
      Particle &p, const MicroXSs &microxs, const Nuclide &nuclide,
      std::vector<std::shared_ptr<NoiseSource>> &noise_sources) const;

  void sample_noise_source_from_fission(
      Particle &p, const std::shared_ptr<Nuclide> nuclide,
      const MicroXSs &microxs,
      std::vector<std::shared_ptr<NoiseSource>> &noise_sources) const;

  void sample_noise_source_from_copy(
      Particle &p,
      std::vector<std::shared_ptr<NoiseSource>> &noise_sources) const;

};  // Transporter

#endif  // MG_TRANSPORTER_H
