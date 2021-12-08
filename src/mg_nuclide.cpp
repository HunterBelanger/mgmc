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
#include <materials/mg_nuclide.hpp>
#include <utils/constants.hpp>
#include <utils/rng.hpp>

MGNuclide::MGNuclide(const std::vector<double> &ebounds,
                     const std::vector<double> &Et,
                     const std::vector<double> &Ea,
                     const std::vector<double> &Ef,
                     const std::vector<double> &nu,
                     const std::vector<double> &chi,
                     const std::vector<std::vector<double>> &Es)
    : energy_bounds(ebounds),
      Et_(Et),
      Ea_(Ea),
      Ef_(Ef),
      nu_prmpt_(nu),
      nu_delyd_(),
      chi_(chi),
      Es_(Es),
      P_delayed_group(),
      delayed_group_decay_constants() {
  // Go through Ef and nu to see if we are fissile or not
  for (std::size_t i = 0; i < Ef_.size(); i++) {
    if (Ef_[i] != 0. && nu_prmpt_[i] != 0.) {
      fissile_ = true;
      break;
    }
  }
}

MGNuclide::MGNuclide(
    const std::vector<double> &ebounds, const std::vector<double> &Et,
    const std::vector<double> &Ea, const std::vector<double> &Ef,
    const std::vector<double> &nu_p, const std::vector<double> &nu_d,
    const std::vector<double> &chi, const std::vector<double> &dgp,
    const std::vector<double> &dgc, const std::vector<std::vector<double>> &Es)
    : energy_bounds(ebounds),
      Et_(Et),
      Ea_(Ea),
      Ef_(Ef),
      nu_prmpt_(nu_p),
      nu_delyd_(nu_d),
      chi_(chi),
      Es_(Es),
      P_delayed_group(dgp),
      delayed_group_decay_constants(dgc) {
  // Go through Ef and nu to see if we are fissile or not
  for (std::size_t i = 0; i < Ef_.size(); i++) {
    if (Ef_[i] != 0. && this->nu_total(0., i) != 0.) {
      fissile_ = true;
      break;
    }
  }
}

bool MGNuclide::fissile() const { return fissile_; }

double MGNuclide::total_xs(double /*E_in*/, std::size_t i) const {
  return Et_[i];
}

double MGNuclide::disappearance_xs(double /*E_in*/, std::size_t i) const {
  return Ea_[i] - Ef_[i];
}

double MGNuclide::fission_xs(double /*E_in*/, std::size_t i) const {
  return Ef_[i];
}

double MGNuclide::nu_total(double /*E_in*/, std::size_t i) const {
  double nu_tot = nu_prmpt_[i];

  if (!nu_delyd_.empty()) nu_tot += nu_delyd_[i];

  return nu_tot;
}

double MGNuclide::nu_prompt(double /*E_in*/, std::size_t i) const {
  return nu_prmpt_[i];
}

double MGNuclide::nu_delayed(double /*E_in*/, std::size_t i) const {
  if (nu_delyd_.empty()) return 0.;

  return nu_delyd_[i];
}

double MGNuclide::reaction_xs(uint32_t /*mt*/, double /*E_in*/,
                              std::size_t /*i*/) const {
  return 0.;
}

double MGNuclide::elastic_xs(double /*E_in*/, std::size_t i) const {
  double Es = 0;

  for (const auto &v : Es_[i]) {
    Es += v;
  }

  return Es;
}

std::size_t MGNuclide::energy_grid_index(double E) const {
  std::size_t i = 0;

  for (i = 0; i < energy_bounds.size() - 1; i++) {
    if (energy_bounds[i] <= E && E < energy_bounds[i + 1]) {
      break;
    }
  }

  return i;
}

std::size_t MGNuclide::num_delayed_groups() const {
  return P_delayed_group.size();
}

double MGNuclide::delayed_group_constant(std::size_t g) const {
  return delayed_group_decay_constants[g];
}

double MGNuclide::delayed_group_probability(std::size_t g, double /*E*/) const {
  return P_delayed_group[g];
}

double MGNuclide::max_energy() const { return energy_bounds.back(); }

double MGNuclide::min_energy() const { return energy_bounds.front(); }

ScatterInfo MGNuclide::sample_scatter(double /*Ein*/, const Direction &u,
                                      std::size_t i, pcg32 &rng) const {
  // Change particle energy
  int ei = RNG::discrete(rng, Es_[i]);
  double E_out = 0.5 * (energy_bounds[ei] + energy_bounds[ei + 1]);

  // Change direction (Isotropic only for now)
  double mu = 2. * RNG::rand(rng) - 1.;
  double phi = 2. * PI * RNG::rand(rng);
  Direction u_out = rotate_direction(u, mu, phi);

  ScatterInfo info;
  info.energy = E_out;
  info.direction = u_out;
  return info;
}

ScatterInfo MGNuclide::sample_prompt_fission(double /*Ein*/, const Direction &u,
                                             std::size_t /*i*/,
                                             pcg32 &rng) const {
  std::function<double()> rngfunc = std::bind(RNG::rand, std::ref(rng));

  // First we sample the the energy index
  int ei = RNG::discrete(rng, chi_);

  // Put fission energy in middle of sampled bin
  double E_out = 0.5 * (energy_bounds[ei] + energy_bounds[ei + 1]);

  // Sample direction from mu, and random phi about z-axis
  double mu = 2. * RNG::rand(rng) - 1.;
  double phi = 2. * PI * RNG::rand(rng);
  Direction u_out = rotate_direction(u, mu, phi);

  ScatterInfo info;
  info.energy = E_out;
  info.direction = u_out;
  return info;
}

ScatterInfo MGNuclide::sample_delayed_fission(double Ein, const Direction &u,
                                              std::size_t /*g*/,
                                              pcg32 &rng) const {
  // Hard code 0 as the variable isn't used in above method for MG nuclides.
  return sample_prompt_fission(Ein, u, 0, rng);
}

BankedParticle MGNuclide::make_fission_neutron(
    Particle &p, std::size_t /*energy_index*/, double P_del,
    std::optional<double> w_noise) const {
  std::function<double()> rngfunc = std::bind(RNG::rand, std::ref(p.rng));

  // First we sample the the energy index
  int ei = RNG::discrete(p.rng, chi_);

  // Put fission energy in middle of sampled bin
  double E_out = 0.5 * (energy_bounds[ei] + energy_bounds[ei + 1]);

  // Sample direction from mu, and random phi about z-axis
  double mu = 2. * RNG::rand(p.rng) - 1.;
  double phi = 2. * PI * RNG::rand(p.rng);
  Direction u = rotate_direction(p.u(), mu, phi);

  // Define weights to be used (will be changed if noise = true).
  double wgt = p.wgt() > 0. ? 1. : -1.;
  double wgt2 = 0.;

  // If we are running in noise mode or sampling the noise source from fission,
  // we need to treat the weights differently.
  if (w_noise) {
    // First, we set wgt and wgt to the be the current weight of the
    // parent particle.
    wgt = p.wgt();
    wgt2 = p.wgt2();

    // Next, we need to see if this is a delayed neutron or not.
    if (rngfunc() < P_del) {
      // We have a delayed neutron. We now need to select a delayed group
      // and get the group decay constant.
      int dgrp = RNG::discrete(p.rng, P_delayed_group);
      double lambda = delayed_group_decay_constants[dgrp];

      // We must now modify this neutron's weight
      std::complex<double> weight{wgt, wgt2};
      double denom = lambda * lambda + w_noise.value() * w_noise.value();
      std::complex<double> mult{lambda / denom,
                                -lambda * w_noise.value() / denom};
      weight *= mult;
      wgt = weight.real();
      wgt2 = weight.imag();
    }
  }

  // Create and return the neutron
  BankedParticle p_new{p.r(),
                       u,
                       E_out,
                       wgt,
                       wgt2,
                       p.history_id(),
                       p.daughter_counter(),
                       p.previous_r(),
                       p.E(),
                       p.Esmp()};

  return p_new;
}

void MGNuclide::scatter(
    Particle &p, std::size_t energy_index, std::optional<double> w_noise,
    std::vector<std::shared_ptr<NoiseSource>> *noise_sources) const {
  // Change particle energy
  int ei = RNG::discrete(p.rng, Es_[energy_index]);
  double E_out = 0.5 * (energy_bounds[ei] + energy_bounds[ei + 1]);

  // Change direction (Isotropic only for now)
  double mu = 2. * RNG::rand(p.rng) - 1.;
  double phi = 2. * PI * RNG::rand(p.rng);
  double ux = std::sqrt(1. - mu * mu) * std::cos(phi);
  double uy = std::sqrt(1. - mu * mu) * std::sin(phi);
  double uz = mu;

  if (noise_sources && w_noise) {
    // Sample new energy and direction for noise particle
    int noise_ei = RNG::discrete(p.rng, Es_[energy_index]);
    double noise_E_out =
        0.5 * (energy_bounds[noise_ei] + energy_bounds[noise_ei + 1]);

    double noise_mu = 2. * RNG::rand(p.rng) - 1.;
    double noise_phi = 2. * PI * RNG::rand(p.rng);
    double noise_ux = std::sqrt(1. - noise_mu * noise_mu) * std::cos(noise_phi);
    double noise_uy = std::sqrt(1. - noise_mu * noise_mu) * std::sin(noise_phi);
    double noise_uz = noise_mu;

    // Make noise particle to bank
    BankedParticle p_noise{p.r(),
                           Direction(noise_ux, noise_uy, noise_uz),
                           noise_E_out,
                           p.wgt(),
                           p.wgt2(),
                           p.history_id(),
                           p.daughter_counter(),
                           p.previous_r(),
                           p.E(),
                           p.Esmp()};

    // Apply weight correction
    std::complex<double> nweight{p.wgt(), p.wgt2()};
    std::complex<double> dEelastic_Eelastic = {0., 0.};
    bool inside_a_noise_source = false;
    for (const auto &ns : *noise_sources) {
      if (ns->is_inside(p.r(), p.u())) {
        inside_a_noise_source = true;
        dEelastic_Eelastic +=
            ns->dEelastic_Eelastic(p.r(), p.u(), p.E(), *w_noise);
      }
    }

    if (inside_a_noise_source) {
      nweight *= dEelastic_Eelastic;

      p_noise.wgt = nweight.real();
      p_noise.wgt2 = nweight.imag();

      p.add_noise_particle(p_noise);
    }
  }

  p.set_energy(E_out);
  p.set_direction(Direction(ux, uy, uz));
}
