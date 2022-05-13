/*=============================================================================*
 * Copyright (C) 2021-2022, Commissariat à l'Energie Atomique et aux Energies
 * Alternatives
 *
 * Contributeur : Hunter Belanger (hunter.belanger@cea.fr)
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
#include <complex>
#include <simulation/transporter.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>

void Transporter::russian_roulette(Particle &p) {
  // Roulette first weight
  if (std::abs(p.wgt()) < settings::wgt_cutoff) {
    double P_kill = 1.0 - (std::abs(p.wgt()) / settings::wgt_survival);
    if (RNG::rand(p.rng) < P_kill)
      p.set_weight(0.);
    else {
      p.set_weight(std::copysign(settings::wgt_survival, p.wgt()));
    }
  }

  // Roulette second weight
  if (std::abs(p.wgt2()) < settings::wgt_cutoff) {
    double P_kill = 1.0 - (std::abs(p.wgt2()) / settings::wgt_survival);
    if (RNG::rand(p.rng) < P_kill)
      p.set_weight2(0.);
    else {
      p.set_weight2(std::copysign(settings::wgt_survival, p.wgt2()));
    }
  }

  // Check if we are still alive
  if (p.wgt() == 0. && p.wgt2() == 0.) p.kill();
}

void Transporter::collision(Particle &p, MaterialHelper &mat,
                            ThreadLocalScores &thread_scores, bool noise,
                            const NoiseMaker *noise_maker) {
  // Score flux collision estimator with Sigma_t
  tallies->score_collision(p, mat, settings::converged);
  double k_col_scr = p.wgt() * mat.vEf(p.E()) / mat.Et(p.E(), noise);
  double mig_dist = (p.r() - p.r_birth()).norm();
  double mig_area_scr =
      p.wgt() * mat.Ea(p.E()) / mat.Et(p.E()) * mig_dist * mig_dist;
  if (!noise) {
    thread_scores.k_col_score += k_col_scr;
    thread_scores.mig_score += mig_area_scr;
  }

  // Sample noise source from the noise maker.
  if (noise_maker) {
    noise_maker->sample_noise_source(p, mat, tallies->keff(),
                                     settings::w_noise);
  }

  // Sample a nuclide for the collision
  std::pair<const std::shared_ptr<Nuclide> &, MicroXSs> nuclide_info =
      mat.sample_nuclide(p.E(), p.rng, noise);
  const std::shared_ptr<Nuclide> &nuclide = nuclide_info.first;
  MicroXSs microxs = nuclide_info.second;

  // Score kabs
  if (!noise) {
    double k_abs_scr =
        p.wgt() * microxs.nu_total * microxs.fission / microxs.total;
    thread_scores.k_abs_score += k_abs_scr;
  }

  // Make all fission neutrons
  make_fission_neutrons(p, microxs, *nuclide, noise);

  // Make a copy of the particle when doing noise transport, modifying
  // the weight of the copy as well.
  if (noise) {
    make_noise_copy(p, microxs);
  }

  // Implicit capture
  p.set_weight(p.wgt() * (1. - (microxs.absorption + microxs.noise_copy) /
                                   microxs.total));
  p.set_weight2(p.wgt2() * (1. - (microxs.absorption + microxs.noise_copy) /
                                     microxs.total));

  // Roulette
  russian_roulette(p);

  // Scatter particle
  if (p.is_alive()) {
    // Perform scatter with nuclide
    nuclide->scatter(p, microxs.energy_index);

    if (p.E() < settings::min_energy) {
      p.kill();
    }
  }
}

void Transporter::make_noise_copy(Particle &p, const MicroXSs &microxs) const {
  std::complex<double> weight_copy{p.wgt(), p.wgt2()};
  if (microxs.noise_copy / microxs.total + RNG::rand(p.rng) >= 1.) {
    std::complex<double> yield{1., -1. / settings::eta};
    weight_copy *= yield;
    p.make_secondary(p.u(), p.E(), weight_copy.real(), weight_copy.imag());
  }
}

void Transporter::make_fission_neutrons(Particle &p, const MicroXSs &microxs,
                                        const Nuclide &nuclide,
                                        bool noise) const {
  double k_abs_scr =
      p.wgt() * microxs.nu_total * microxs.fission / microxs.total;

  // Determine the number of neutrons to produce
  int n_new = 0;
  switch (settings::mode) {
    case settings::SimulationMode::K_EIGENVALUE:
      // In k-eivenvalue simulations, we normalize particle production
      // by keff so that it stays more or less constant.
      n_new =
          std::floor(std::abs(k_abs_scr) / tallies->kcol() + RNG::rand(p.rng));
      break;

    case settings::SimulationMode::FIXED_SOURCE:
      // In fixed-source problems, we don't normalize the number of
      // particles. Problem MUST BE SUB CRITICAL !
      n_new = std::floor(std::abs(k_abs_scr) + RNG::rand(p.rng));
      break;

    case settings::SimulationMode::NOISE:
      if (!noise) {
        // This is a power iteration generation. Run as with K_EIGENVALUE
        n_new = std::floor(std::abs(k_abs_scr) / tallies->kcol() +
                           RNG::rand(p.rng));
      } else {
        // When transporting noise particles, we don't normalize the number of
        // generated particles by the weight ! We also scale by kcol from the
        // previous power-iteration generation, to make the problem "critical".
        n_new = std::floor((microxs.nu_total * microxs.fission /
                            (microxs.total * tallies->keff())) +
                           RNG::rand(p.rng));
      }
      break;
  }

  // Probability of a fission neutron being a delayed neutron.
  double P_delayed = microxs.nu_delayed / microxs.nu_total;

  // If we are are not looking at a noise particle, then the noise frequency
  // will not be given, which tells the make_fission_neutron method to not use
  // the parent weights, and don't consider the delayed noise factor. Otherwise,
  // w_noise is the sought noise frequency, and this will indicate that noise
  // particle magic should be done when making the fission neutron.
  std::optional<double> w_noise = std::nullopt;
  if (noise) w_noise = settings::w_noise;

  for (int i = 0; i < n_new; i++) {
    auto fiss_particle = nuclide.make_fission_neutron(p, microxs.energy_index,
                                                      P_delayed, w_noise);

    switch (settings::mode) {
      case settings::SimulationMode::FIXED_SOURCE:
        p.make_secondary(fiss_particle.u, fiss_particle.E, fiss_particle.wgt,
                         fiss_particle.wgt2);
        break;

      case settings::SimulationMode::K_EIGENVALUE:
        p.add_fission_particle(fiss_particle);
        break;

      case settings::SimulationMode::NOISE:
        if (!noise || settings::inner_generations) {
          p.add_fission_particle(fiss_particle);
        } else {
          p.make_secondary(fiss_particle.u, fiss_particle.E, fiss_particle.wgt,
                           fiss_particle.wgt2);
        }
        break;
    }
  }
}
