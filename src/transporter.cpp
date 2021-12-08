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

void Transporter::collision(
    Particle &p, MaterialHelper &mat, ThreadLocalScores &thread_scores,
    bool noise, std::vector<std::shared_ptr<NoiseSource>> *noise_sources) {
  // Score flux collision estimator with Sigma_t
  tallies->score_collision(p, mat, settings::converged);
  double k_col_scr = p.wgt() * mat.vEf(p.E()) / mat.Et(p.E(), noise);
  if (!noise) {
    thread_scores.k_col_score += k_col_scr;
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

  if (noise_sources && settings::branchless_noise) {
    sample_branchless_noise(p, microxs, *nuclide, *noise_sources);
  } else if (noise_sources && !settings::branchless_noise) {
    // Sample fission portion of noise source.
    sample_noise_source_from_fission(p, nuclide, microxs, *noise_sources);
    sample_noise_source_from_copy(p, *noise_sources);
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
    // Noise frequency
    std::optional<double> w_noise = std::nullopt;
    if (noise || noise_sources) w_noise = settings::w_noise;
    if (settings::branchless_noise) {
      w_noise = std::nullopt;
      noise_sources = nullptr;
    }

    // Perform scatter with nuclide
    nuclide->scatter(p, microxs.energy_index, w_noise, noise_sources);

    if (p.E() < settings::min_energy) {
      p.kill();
    }
  }
}

void Transporter::sample_branchless_noise(
    Particle &p, const MicroXSs &microxs, const Nuclide &nuclide,
    std::vector<std::shared_ptr<NoiseSource>> &noise_sources) const {
  // When doing the noise sampling, we have 3+Ng terms, where Ng is the number
  // of delayed fission groups.
  const std::size_t Ng{nuclide.num_delayed_groups()};
  std::size_t Nt = 3 + Ng;

  // Initialize vector to contain the info
  std::vector<std::complex<double>> yield_sigma(Nt, {0., 0.});
  std::vector<double> mod(Nt, 0.);
  double sum_mods = 0.;

  // Get xs ratios
  std::complex<double> dEt_Et = {0., 0.};
  std::complex<double> dEs_Es = {0., 0.};
  std::complex<double> dEf_Ef = {0., 0.};
  bool inside_a_noise_source = false;
  for (const auto &ns : noise_sources) {
    if (ns->is_inside(p.r(), p.u())) {
      inside_a_noise_source = true;
      dEt_Et += ns->dEt_Et(p.r(), p.u(), p.E(), settings::w_noise);
      dEs_Es += ns->dEelastic_Eelastic(p.r(), p.u(), p.E(), settings::w_noise);
      dEf_Ef += ns->dEf_Ef(p.r(), p.u(), p.E(), settings::w_noise);
    }
  }

  // Not in noise source, don't make noise particles
  if (!inside_a_noise_source) return;

  // Copy Info
  yield_sigma[0] = -dEt_Et * microxs.total;

  // Scatter Info
  yield_sigma[1] = dEs_Es * (microxs.total - microxs.absorption);

  // Promp Fission Info
  double vt = microxs.nu_total;
  double vd = microxs.nu_delayed;
  double vp = vt - vd;
  yield_sigma[2] = vp * dEf_Ef * microxs.fission;

  // Delayed Fission Info
  double w = settings::w_noise;
  for (std::size_t g = 0; g < Ng; g++) {
    double Pg = nuclide.delayed_group_probability(g, p.E());
    double lambda = nuclide.delayed_group_constant(g);
    double denom = lambda * lambda + w * w;
    yield_sigma[3 + g] = {lambda * lambda / denom, -lambda * w / denom};
    yield_sigma[3 + g] *= vd * dEf_Ef * microxs.fission * Pg;
  }

  // Get all Mods
  for (std::size_t i = 0; i < Nt; i++) {
    mod[i] = std::abs(yield_sigma[i]);
    sum_mods += mod[i];
  }

  // Sample reaction
  const int j{RNG::discrete(p.rng, mod)};
  double m = sum_mods / microxs.total;
  std::complex<double> m_j = m * yield_sigma[j] / mod[j];
  double E_out = p.E();
  Direction u_out = p.u();

  if (j >= static_cast<int>(mod.size())) {
    // Shouldn't get here ! RNG::discrete is rather safe as I use
    // std::discrete_distribution, but because I don't have a default
    // check in the swithc, it's here anyway.
    std::string mssg = "Invalid index when sampling branchless noise source.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  switch (j) {
    case 0:
      // Copy
      // Nothing to do, E_out and u_out initialized as the energy and
      // direction from p already.
      break;

    case 1: {
      // Scatter
      auto info =
          nuclide.sample_scatter(p.E(), p.u(), microxs.energy_index, p.rng);
      E_out = info.energy;
      u_out = info.direction;
      m_j *= info.yield;
      break;
    }

    case 2: {
      // Prompt Fission
      auto info = nuclide.sample_prompt_fission(p.E(), p.u(),
                                                microxs.energy_index, p.rng);
      E_out = info.energy;
      u_out = info.direction;
      break;
    }

    default: {
      // Delayed Fission
      auto info = nuclide.sample_delayed_fission(p.E(), p.u(), j, p.rng);
      E_out = info.energy;
      u_out = info.direction;
      break;
    }
  }

  std::complex<double> wgt = {p.wgt(), p.wgt2()};
  wgt *= m_j;

  // Create and return the neutron
  BankedParticle p_noise{p.r(),
                         u_out,
                         E_out,
                         wgt.real(),
                         wgt.imag(),
                         p.history_id(),
                         p.daughter_counter(),
                         p.previous_r(),
                         p.E(),
                         p.Esmp()};

  p.add_noise_particle(p_noise);
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

void Transporter::sample_noise_source_from_copy(
    Particle &p,
    std::vector<std::shared_ptr<NoiseSource>> &noise_sources) const {
  std::complex<double> weight_copy{p.wgt(), p.wgt2()};

  // Get the weight modifier
  std::complex<double> dEt_Et = {0., 0.};
  bool inside_a_noise_source = false;
  for (const auto &ns : noise_sources) {
    if (ns->is_inside(p.r(), p.u())) {
      inside_a_noise_source = true;
      dEt_Et += ns->dEt_Et(p.r(), p.u(), p.E(), settings::w_noise);
    }
  }

  if (inside_a_noise_source == false) {
    // Not inside any source. Don't make noise
    return;
  }

  // Negative needed for source sampling
  weight_copy *= -dEt_Et;

  // Create and return the neutron
  BankedParticle p_noise{p.r(),
                         p.u(),
                         p.E(),
                         weight_copy.real(),
                         weight_copy.imag(),
                         p.history_id(),
                         p.daughter_counter(),
                         p.previous_r(),
                         p.E(),
                         p.Esmp()};

  p.add_noise_particle(p_noise);
}

void Transporter::sample_noise_source_from_fission(
    Particle &p, const std::shared_ptr<Nuclide> nuclide,
    const MicroXSs &microxs,
    std::vector<std::shared_ptr<NoiseSource>> &noise_sources) const {
  int n_new = 0;

  // We don't multiply by the weight when getting the number of noise
  // particles to make, as in the make_fission_neutron method, if we
  // pass the noise frequency, it will automatically use the incident
  // particle weight in the daughter particle weight.
  double k_abs_scr = microxs.nu_total * microxs.fission / microxs.total;
  n_new = std::floor(k_abs_scr / tallies->keff() + RNG::rand(p.rng));

  // Probability of a fission neutron being a delayed neutron.
  double P_delayed = microxs.nu_delayed / microxs.nu_total;

  // Noise frequency
  double w_noise = settings::w_noise;

  // Get the weight modifier
  std::complex<double> dEf_Ef = {0., 0.};
  bool inside_a_noise_source = false;
  for (const auto &ns : noise_sources) {
    if (ns->is_inside(p.r(), p.u())) {
      inside_a_noise_source = true;
      dEf_Ef += ns->dEf_Ef(p.r(), p.u(), p.E(), w_noise);
    }
  }

  if (inside_a_noise_source == false) {
    // Don't make any particles !
    return;
  }

  for (int i = 0; i < n_new; i++) {
    auto fiss_noise_particle = nuclide->make_fission_neutron(
        p, microxs.energy_index, P_delayed, w_noise);

    // Modify weight
    std::complex<double> fnp_weight{fiss_noise_particle.wgt,
                                    fiss_noise_particle.wgt2};
    fnp_weight *= dEf_Ef;
    fiss_noise_particle.wgt = fnp_weight.real();
    fiss_noise_particle.wgt2 = fnp_weight.imag();

    // Save BankedParticle
    p.add_noise_particle(fiss_noise_particle);
  }
}
