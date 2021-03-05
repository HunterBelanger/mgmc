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
#include <omp.h>

#include <simulation/carter_tracker.hpp>
#include <simulation/tracker.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>

std::vector<Particle> CarterTracker::transport(
    std::vector<Particle>& bank, std::vector<std::shared_ptr<RNG>>& rngs) {
  // Vector to contain all fission daughters for all threads
  std::vector<Particle> fission_neutrons;
  std::vector<Particle> Bank = bank;
  int n_particles = static_cast<int>(Bank.size());

  while (n_particles > 0) {
#pragma omp parallel
    {
      // Get random generator for thread
      int thrd = 0;
#ifdef _OPENMP
      thrd = omp_get_thread_num();
#endif
      std::shared_ptr<RNG> rng = rngs[thrd];

      // Thread local storage
      std::vector<Particle> thread_fissions;
      thread_fissions.reserve(settings->nparticles / rngs.size());
      std::vector<Particle> thread_splits;

// int k;

// Transport all particles in for thread
#pragma omp for schedule(static)
      for (size_t n = 0; n < Bank.size(); n++) {
        // Particle and its personal tracker
        Particle p = Bank[n];
        Tracker trkr(p.r(), p.u());

        // If we got lost, kill the particle
        if (trkr.is_lost()) {
          std::stringstream mssg;
          mssg << "Particle become lost at " << p.r() << ".";
          warning(mssg.str(), __FILE__, __LINE__);
          p.kill();
        }

        while (p.is_alive()) {
          bool virtual_collision = true;
          while (virtual_collision) {
            double d_coll = rng->exponential(Esmp[p.E()]);
            geometry::Boundary surf_bound = trkr.boundary();

            // std::cout << p.r() << "  " << p.u() << "\n";
            // std::cout << "d_col = "<<d_coll<<", d_surf =
            // "<<surf_bound.distance()<<"\n";

            if (d_coll < surf_bound.distance) {
              p.move(d_coll);
              trkr.move(d_coll);
              trkr.get_current();
              double Et = trkr.material()->Et(p.r(), p.E());
              if (Esmp[p.E()] >= Et) {
                if (rng->rand() < Et / Esmp[p.E()]) virtual_collision = false;
              } else {
                double D = Et / (2. * Et - Esmp[p.E()]);
                double F = Et / (D * Esmp[p.E()]);

                if ((D - rng->rand()) > 0.) {
                  p.set_weight(p.wgt() * F);
                  virtual_collision = false;
                } else {
                  p.set_weight(-p.wgt() * F);
                }
              }
            } else {
              if (surf_bound.boundary_type == BoundaryType::Vacuum) {
                p.kill();
                tallies->score_leak(p.wgt());
                virtual_collision = false;
                // std::cout << " Particle has leaked...\n";
              } else if (surf_bound.boundary_type == BoundaryType::Reflective) {
                // geometry::do_reflection(p, surf_bound);
                // trkr.set_r(p.r()); trkr.set_u(p.u());
                // trkr.restart_get_current();
                trkr.do_reflection(p, surf_bound);
              }  // Reflection
            }    // Boundary interaction
          }      // While virtual collision

          // If we got lost, kill the particle
          if (trkr.is_lost()) {
            std::stringstream mssg;
            mssg << "Particle become lost at " << p.r() << ".";
            warning(mssg.str(), __FILE__, __LINE__);
            p.kill();
          }

          if (p.is_alive()) {  // real collision
            std::shared_ptr<Material> mat = trkr.material();
            double Et = mat->Et(p.r(), p.E());
            double Ea = mat->Ea(p.r(), p.E());
            double Ef = mat->Ef(p.r(), p.E());
            double nu = mat->nu(p.r(), p.E());

            // Do scores
            // std::cout << " New location: " << p.r() << "\n";
            tallies->score_collision(p, mat, settings->converged);

            // std::cout<<" p.wgt = "<<p.wgt()<<", p.E = "<<p.E()<<", k_col_scr
            // = "<<k_col_scr << "\n";

            // Get new fission neutrons
            double k_col_scr = p.wgt() * nu * Ef / Et;
            int n_new =
                std::floor(std::abs(k_col_scr) / tallies->keff() + rng->rand());
            // std::cout << " n_new = " << n_new << "\n";
            for (int i = 0; i < n_new; i++) {
              thread_fissions.push_back(fission_neutron(p, mat, rng));
            }

            // Implicit capture
            p.set_weight(p.wgt() * (1. - (Ea / Et)));
            // std::cout << " New p.wgt = " << p.wgt() << "\n";

            // Scatter particle
            scatter_particle(p, mat, rng);
            trkr.set_u(p.u());

            // Roulette
            russian_roulette(p, rng);
            // if(p.is_alive() == false) std::cout << " particle killed by
            // roulette...\n"; std::cout << "\n";

          }  // If alive for real collision

          // Split partile if needed
          if (p.is_alive() and (std::abs(p.wgt()) >= settings->wgt_split)) {
            double n_new = std::ceil(std::abs(p.wgt()));
            p.set_weight(p.wgt() / n_new);
            for (int np = 0; np < static_cast<int>(n_new - 1); np++) {
              Particle p_daughter(p.r(), p.u(), p.E(), p.wgt());
              thread_splits.push_back(p_daughter);
            }
          }

        }  // While alive
        // std::cout << " END PARTICLE HISTORY\n\n";
        // std::cin >> k;
      }  // For all particles

      // Save random number generator
      rngs[thrd] = rng;

#pragma omp single
      { Bank.clear(); }

#pragma omp barrier
#pragma omp critical
      {
        fission_neutrons.insert(std::end(fission_neutrons),
                                std::begin(thread_fissions),
                                std::end(thread_fissions));

        Bank.insert(std::end(Bank), std::begin(thread_splits),
                    std::end(thread_splits));
      }

    }  // Parallel

    n_particles = static_cast<int>(Bank.size());

  }  // While still particles to transport

  return fission_neutrons;
}

void CarterTracker::russian_roulette(Particle& p, std::shared_ptr<RNG>& rng) {
  if (std::abs(p.wgt()) < settings->wgt_cutoff) {
    double P_kill = 1.0 - (std::abs(p.wgt()) / settings->wgt_survival);
    if (rng->rand() < P_kill)
      p.kill();
    else {
      if (p.wgt() > 0.)
        p.set_weight(settings->wgt_survival);
      else
        p.set_weight(-settings->wgt_survival);
    }
  }
}

void CarterTracker::scatter_particle(Particle& p, std::shared_ptr<Material> mat,
                                     std::shared_ptr<RNG>& rng) {
  // Change particle energy
  p.set_energy(rng->discrete(mat->Es(p.r(), p.E())));

  // Change direction (Isotropic only for now)
  double mu = 2. * rng->rand() - 1.;
  double phi = 2. * PI * rng->rand();
  double ux = std::sqrt(1. - mu * mu) * std::cos(phi);
  double uy = std::sqrt(1. - mu * mu) * std::sin(phi);
  double uz = mu;

  p.set_direction(Direction(ux, uy, uz));
}

Particle CarterTracker::fission_neutron(Particle& p,
                                        std::shared_ptr<Material> mat,
                                        std::shared_ptr<RNG>& rng) {
  // Ge weight
  double w;
  if (p.wgt() > 0.)
    w = 1.;
  else
    w = -1.;

  // Get energy
  int E = rng->discrete(mat->chi(p.r(), p.E()));

  // Get direction
  double mu = 2. * rng->rand() - 1.;
  double phi = 2. * PI * rng->rand();
  double ux = std::sqrt(1. - mu * mu) * std::cos(phi);
  double uy = std::sqrt(1. - mu * mu) * std::sin(phi);
  double uz = mu;
  Direction u(ux, uy, uz);

  Particle p_new(p.r(), u, E, w);
  p_new.parents_previous_energy = p.E();
  p_new.parents_previous_position = p.previous_r();
  p_new.Esmp_parent = Esmp[p.E()];

  return p_new;
}
