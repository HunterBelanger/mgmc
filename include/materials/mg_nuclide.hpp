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
#ifndef MG_NUCLIDE_H
#define MG_NUCLIDE_H

#include <materials/nuclide.hpp>
#include <vector>

class MGNuclide : public Nuclide {
 public:
  MGNuclide(const std::vector<double> &ebounds, const std::vector<double> &Et,
            const std::vector<double> &Ea, const std::vector<double> &Ef,
            const std::vector<double> &nu, const std::vector<double> &chi,
            const std::vector<std::vector<double>> &Es);

  MGNuclide(const std::vector<double> &ebounds, const std::vector<double> &Et,
            const std::vector<double> &Ea, const std::vector<double> &Ef,
            const std::vector<double> &nu_p, const std::vector<double> &nu_d,
            const std::vector<double> &chi, const std::vector<double> &dgp,
            const std::vector<double> &dgc,
            const std::vector<std::vector<double>> &Es);

  bool fissile() const override final;
  double total_xs(double E_in, std::size_t i) const override final;
  double disappearance_xs(double E_in, std::size_t i) const override final;
  double fission_xs(double E_in, std::size_t i) const override final;
  double nu_total(double E_in, std::size_t i) const override final;
  double nu_prompt(double E_in, std::size_t i) const override final;
  double nu_delayed(double E_in, std::size_t i) const override final;
  double reaction_xs(uint32_t mt, double E_in, size_t i) const override final;
  double elastic_xs(double E_in, std::size_t i) const override final;
  std::size_t energy_grid_index(double E) const override final;
  std::size_t num_delayed_groups() const override final;
  double delayed_group_constant(std::size_t g) const override final;
  double delayed_group_probability(std::size_t g,
                                   double E) const override final;
  ScatterInfo sample_scatter(double Ein, const Direction &u, std::size_t i,
                             pcg32 &rng) const override final;
  ScatterInfo sample_prompt_fission(double Ein, const Direction &u,
                                    std::size_t i,
                                    pcg32 &rng) const override final;
  ScatterInfo sample_delayed_fission(double Ein, const Direction &u,
                                     std::size_t g,
                                     pcg32 &rng) const override final;

  double max_energy() const override final;
  double min_energy() const override final;

  BankedParticle make_fission_neutron(
      Particle &p, std::size_t energy_index, double P_del,
      std::optional<double> w_noise = std::nullopt) const override final;
  void scatter(Particle &p, std::size_t energy_index,
               std::optional<double> w_noise = std::nullopt,
               std::vector<std::shared_ptr<NoiseSource>> *noise_sources =
                   nullptr) const override final;

 private:
  std::vector<double> energy_bounds;
  std::vector<double> Et_;
  std::vector<double> Ea_;
  std::vector<double> Ef_;
  std::vector<double> nu_prmpt_;
  std::vector<double> nu_delyd_;
  std::vector<double> chi_;
  std::vector<std::vector<double>> Es_;
  std::vector<double> P_delayed_group;
  std::vector<double> delayed_group_decay_constants;
  bool fissile_ = false;
};

#endif
