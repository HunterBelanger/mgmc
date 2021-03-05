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
#ifndef MG_SIMULATION_H
#define MG_SIMULATION_H

#include <geometry/geometry.hpp>
#include <simulation/entropy.hpp>
#include <simulation/source.hpp>
#include <simulation/tallies.hpp>
#include <simulation/transporter.hpp>
#include <sstream>
#include <utils/rng.hpp>
#include <utils/settings.hpp>

class Simulation {
 public:
  Simulation(std::shared_ptr<Tallies> i_t, std::shared_ptr<Settings> i_s,
             std::shared_ptr<Transporter> i_tr,
             std::vector<std::shared_ptr<Source>> srcs)
      : tallies{i_t},
        settings{i_s},
        transporter{i_tr},
        sources{srcs},
        rngs{},
        p_pre_entropy(nullptr),
        n_pre_entropy(nullptr),
        t_pre_entropy(nullptr),
        p_post_entropy(nullptr),
        n_post_entropy(nullptr),
        t_post_entropy(nullptr){};
  virtual ~Simulation() = default;

  virtual void initialize() = 0;
  virtual void run() = 0;

  virtual void premature_kill() = 0;

  // Methods to set entropies
  void set_p_pre_entropy(std::shared_ptr<Entropy> entrpy) {
    p_pre_entropy = entrpy;
  }
  void set_n_pre_entropy(std::shared_ptr<Entropy> entrpy) {
    n_pre_entropy = entrpy;
  }
  void set_t_pre_entropy(std::shared_ptr<Entropy> entrpy) {
    t_pre_entropy = entrpy;
  }

  void set_p_post_entropy(std::shared_ptr<Entropy> entrpy) {
    p_post_entropy = entrpy;
  }
  void set_n_post_entropy(std::shared_ptr<Entropy> entrpy) {
    n_post_entropy = entrpy;
  }
  void set_t_post_entropy(std::shared_ptr<Entropy> entrpy) {
    t_post_entropy = entrpy;
  }

  bool signaled = false;

 protected:
  std::shared_ptr<Tallies> tallies;
  std::shared_ptr<Settings> settings;
  std::shared_ptr<Transporter> transporter;
  std::vector<std::shared_ptr<Source>> sources;
  std::vector<std::shared_ptr<RNG>> rngs;

  // All entropy bins possible
  std::shared_ptr<Entropy> p_pre_entropy = nullptr;
  std::shared_ptr<Entropy> n_pre_entropy = nullptr;
  std::shared_ptr<Entropy> t_pre_entropy = nullptr;

  std::shared_ptr<Entropy> p_post_entropy = nullptr;
  std::shared_ptr<Entropy> n_post_entropy = nullptr;
  std::shared_ptr<Entropy> t_post_entropy = nullptr;

};  // Simulation

#endif  // MG_SIMULATION_H
