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
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <PapillonNDL/cross_section.hpp>
#include <algorithm>
#include <iomanip>
#include <materials/material.hpp>
#include <materials/material_helper.hpp>
#include <simulation/delta_tracker.hpp>
#include <simulation/tracker.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>
#include <vector>

DeltaTracker::DeltaTracker(std::shared_ptr<Tallies> i_t)
    : Transporter(i_t), EGrid(nullptr), Emaj(nullptr) {
  Output::instance()->write(" Finding majorant cross sections.\n");
  // Must first create a unionized energy grid. How this is done depends on
  // whether or not we are in continuous energy or multi-group mode.
  if (settings::energy_mode == settings::EnergyMode::CE) {
    std::string mssg = "Continuous-Energy mode node supported.";
    fatal_error(mssg, __FILE__, __LINE__);
  } else {
    // We are in multi-group mode. Here, the energy-bounds are kept in the
    // settings, so we can construct something with that

    std::vector<double> egrid;
    egrid.push_back(settings::energy_bounds[0]);
    if (egrid.front() == 0.) egrid.front() = 1.E-11;

    for (size_t i = 1; i < settings::energy_bounds.size() - 1; i++) {
      egrid.push_back(settings::energy_bounds[i]);
      egrid.push_back(settings::energy_bounds[i]);
    }

    egrid.push_back(settings::energy_bounds.back());

    // This now has created the vector egrid which will look something like this
    // for the case of 5 energy groups.
    // [0., 1.,   1., 2.,   2., 3.,   3., 4.,   4., 5.]
    // This works, because the energy of multi-group particles should always be
    // inbetween the bounds for the group.

    // Now we need to make a vector which will contian the majorant cross
    // cross sections for each group.
    std::vector<double> maj_xs(egrid.size(), 0.);

    // We loop through materials
    for (const auto &material : materials) {
      // Then we loop through energies
      MaterialHelper mat(material.second, 1.);

      for (uint32_t g = 0; g < settings::ngroups; g++) {
        // Get the energy at the mid-point for the group
        size_t i = g * 2;
        double Eg = 0.5 * (egrid[i] + egrid[i + 1]);

        double xs = mat.Et(Eg);
        if (xs > maj_xs[i]) {
          maj_xs[i] = xs;
          maj_xs[i + 1] = xs;
        }
      }
    }

    // Now we construct the energy grid and majorant
    EGrid = std::make_shared<pndl::EnergyGrid>(egrid, settings::ngroups);
    Emaj = std::make_shared<pndl::CrossSection>(maj_xs, EGrid, 0);
  }
}

std::vector<BankedParticle> DeltaTracker::transport(
    std::vector<Particle> &bank, bool noise,
    std::vector<BankedParticle> *noise_bank,
    std::vector<std::shared_ptr<NoiseSource>> *noise_sources) {
#ifdef _OPENMP
#pragma omp parallel
#endif
  {
    // Thread local storage
    ThreadLocalScores thread_scores;

// Transport all particles in for thread
#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (size_t n = 0; n < bank.size(); n++) {
      // Particle and its personal tracker
      Particle &p = bank[n];
      Tracker trkr(p.r(), p.u());

      // If we got lost, kill the particle
      if (trkr.is_lost()) {
        std::stringstream mssg;
        mssg << "Particle become lost at " << p.r() << ", ";
        mssg << " u = " << p.u() << ", token = " << trkr.surface_token();
        warning(mssg.str(), __FILE__, __LINE__);
        p.kill();
      }
      // Only make helper if we aren't lost, to make sure that material isn't
      // a nullptr
      MaterialHelper mat(trkr.material(), p.E());

      // auto bound = trkr.boundary();
      while (p.is_alive()) {
        bool had_collision = false;
        auto maj_indx = EGrid->get_lower_index(p.E());
        double Emajorant =
            Emaj->evaluate(p.E(), maj_indx) + mat.Ew(p.E(), noise);
        p.set_Esmp(Emajorant);  // Sampling XS saved for cancellation
        double d_coll = RNG::exponential(p.rng, Emajorant);
        auto bound = trkr.boundary();

        if (bound.distance < d_coll) {
          if (bound.boundary_type == BoundaryType::Vacuum) {
            p.kill();
            thread_scores.leakage_score += p.wgt();
            Position r_leak = p.r() + bound.distance * p.u();
          } else if (bound.boundary_type == BoundaryType::Reflective) {
            trkr.do_reflection(p, bound);
          } else {
            fatal_error("Help me, how did I get here ?", __FILE__, __LINE__);
          }
        } else {
          // Update Position
          p.move(d_coll);
          trkr.move(d_coll);
          trkr.get_current();
          mat.set_material(trkr.material(), p.E());

          // Get true cross section here
          double Et = mat.Et(p.E(), noise);

          if (Et - Emajorant > 1.E-10) {
#ifdef _OPENMP
#pragma omp critical
#endif
            {
              std::cout << " E = " << std::setprecision(18) << p.E() << "\n";
              std::cout << " Et   = " << std::setprecision(18) << Et << "\n";
              std::cout << " Emaj = " << std::setprecision(18) << Emajorant
                        << "\n";
              std::string mssg = "Total cross section excedeed majorant.";
              fatal_error(mssg, __FILE__, __LINE__);
            }
          }

          if (RNG::rand(p.rng) < (Et / Emajorant)) {
            // Flag real collision
            had_collision = true;
          }
        }

        if (p.is_alive() && had_collision) {  // real collision
          collision(p, mat, thread_scores, noise, noise_sources);
          trkr.set_u(p.u());
        }  // If alive for real collision

        if (!p.is_alive()) {
          // Attempt a resurection
          p.resurect();

          if (p.is_alive()) {
            trkr.set_r(p.r());
            trkr.set_u(p.u());
            trkr.restart_get_current();
            mat.set_material(trkr.material(), p.E());
          } else if (settings::rng_stride_warnings) {
            // History is truly dead.
            // Check if we went past the particle stride.
            uint64_t n_rng_calls = p.number_of_rng_calls();
            if (n_rng_calls > settings::rng_stride) {
              // This isn't really a good thing. We should
              // write a warning.
              std::string mssg = "History " + std::to_string(p.history_id()) +
                                 " overran the RNG stride.";
              warning(mssg, __FILE__, __LINE__);
            }
          }
        }
      }  // While alive
    }    // For all particles

    // Send all thread local scores to tallies instance
    tallies->score_k_col(thread_scores.k_col_score);
    tallies->score_k_abs(thread_scores.k_abs_score);
    tallies->score_k_trk(thread_scores.k_trk_score);
    tallies->score_k_tot(thread_scores.k_tot_score);
    tallies->score_leak(thread_scores.leakage_score);
    thread_scores.k_col_score = 0.;
    thread_scores.k_abs_score = 0.;
    thread_scores.k_trk_score = 0.;
    thread_scores.k_tot_score = 0.;
    thread_scores.leakage_score = 0.;
  }  // Parallel

  // Vector to contain all fission daughters for all threads
  std::vector<BankedParticle> fission_neutrons;

  // Empty all particle fission banks into the main one
  for (auto &p : bank) {
    p.empty_fission_bank(fission_neutrons);
  }

  if (noise_bank && noise_sources) {
    for (auto &p : bank) {
      p.empty_noise_bank(*noise_bank);
    }
  }

  // Can now clear the old bank
  bank.clear();

  return fission_neutrons;
}
