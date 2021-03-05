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

#include <fstream>
#include <iomanip>
#include <iostream>
#include <simulation/power_iterator.hpp>
#include <utils/mt19937.hpp>
#include <utils/output.hpp>
#include <utils/pcg.hpp>

void PowerIterator::initialize() {
  std::shared_ptr<Output> out = Output::instance();
  out->write(" Seeding random number generators...\n");
  // Initialize random number generators
  int nthrds = 1;
#ifdef _OPENMP
  nthrds = omp_get_max_threads();
#endif
  for (int i = 0; i < nthrds; i++) {
    // Default use PCG for RNG
    if (settings->use_pcg)
      rngs.push_back(std::make_shared<PCG>(settings->rng_seed + i + 5));
    // Otherwise use Mersenne Twister
    else
      rngs.push_back(std::make_shared<MT19937>(settings->rng_seed + i + 5));
  }

  // If not using source file
  if (settings->load_source_file == false) {
    // Vector of source weights
    std::vector<double> wgts;
    for (size_t i = 0; i < sources.size(); i++)
      wgts.push_back(sources[i]->wgt());

    // Generate source particles
    out->write(" Generating source particles...\n");
    std::vector<Particle> source_particles;
    std::shared_ptr<RNG> rng = rngs[0];
    for (int i = 0; i < settings->nparticles; i++) {
      int indx = rng->discrete(wgts);
      source_particles.push_back(sources[indx]->generate_particle(rng));
    }
    bank = source_particles;
  } else {
    // Load source file
    NDArray<double> source =
        NDArray<double>::load(settings->in_source_file_name);
    // Get number of particles
    size_t Nprt = source.shape()[0];

    std::shared_ptr<RNG> rng = rngs[0];
    double tot_wgt = 0.0;
    for (size_t i = 0; i < Nprt; i++) {
      double x = source[i * 5 + 0];
      double y = source[i * 5 + 1];
      double z = source[i * 5 + 2];
      int E = static_cast<int>(source[i * 5 + 3]);
      double w = source[i * 5 + 4];

      double mu = 2. * rng->rand() - 1.;
      double phi = 2. * PI * rng->rand();
      double ux = std::sqrt(1. - mu * mu) * std::cos(phi);
      double uy = std::sqrt(1. - mu * mu) * std::sin(phi);
      double uz = mu;

      bank.push_back({{x, y, z}, {ux, uy, uz}, E, w});
      tot_wgt += w;
    }

    out->write(" Total Weight of System: " + std::to_string(tot_wgt) + "\n");
    settings->nparticles = std::round(tot_wgt);
    tallies->set_total_weight(tot_wgt);
  }
}

void PowerIterator::print_header() {
  std::shared_ptr<Output> out = Output::instance();

  // Change the header for the output based on wether or not the entropy is
  // being calculated
  out->write("\n");
  if (t_pre_entropy) {
    out->write(
        "                                              Pre-Cancelation "
        "Entropies             Post-Cancelation Entropies\n");
    out->write(
        "                                        "
        "------------------------------------   "
        "------------------------------------\n");
    out->write(
        "  Gen      k         kavg +/- err         Positive     Nevative       "
        " Total     Positive     Negative        Total      Nnet      Ntot     "
        " Npos      Nneg      Wnet      Wtot      Wpos      Wneg\n");
    out->write(
        " ---------------------------------------------------------------------"
        "----------------------------------------------------------------------"
        "-------------------------------------------------------\n");
    //             1120   1.23456   1.23456 +/-
    //             0.00023   3.283E+002   3.283E+002   13.283E+002   3.283E+002
    //             3.283E+002   3.283E+002   234567   1234567   1234567 1234567
  } else {
    out->write(
        "  Gen      k         kavg +/- err          Nnet      Ntot      Npos   "
        "   Nneg      Wnet      Wtot      Wpos      Wneg\n");
    out->write(
        " ---------------------------------------------------------------------"
        "-----------------------------------------------\n");
    //             1120   1.23456   1.23456 +/- 0.00023   1234567   1234567
    //             1234567   1234567
  }
}

void PowerIterator::run() {
  std::shared_ptr<Output> out = Output::instance();
  out->write(" Running k Eigenvalue Problem...\n");
  out->write(" NPARTICLES: " + std::to_string(settings->nparticles) + ", ");
  out->write(" NGENERATIONS: " + std::to_string(settings->ngenerations) + ",");
  out->write(" NIGNORED: " + std::to_string(settings->nignored) + "\n");

  // Zero all entropy bins
  if (t_pre_entropy) {
    p_pre_entropy->zero();
    n_pre_entropy->zero();
    t_pre_entropy->zero();

    p_post_entropy->zero();
    n_post_entropy->zero();
    t_post_entropy->zero();
  }

  for (int g = 1; g <= settings->ngenerations; g++) {
    std::vector<Particle> next_gen = transporter->transport(bank, rngs);

    // Get new keff
    tallies->calc_gen_values();

    // Store gen if passed ignored
    if (settings->converged) tallies->record_generation();

    // Zero tallies for next generation
    tallies->clear_generation(settings->converged);

    //---------------------------------------------------------------
    // Do all Pre-Cancelation entropy calculations
    if (t_pre_entropy) {
#pragma omp parallel for schedule(static)
      for (size_t i = 0; i < next_gen.size(); i++) {
        p_pre_entropy->add_point(next_gen[i].r(), next_gen[i].wgt());
        n_pre_entropy->add_point(next_gen[i].r(), next_gen[i].wgt());
        t_pre_entropy->add_point(next_gen[i].r(), next_gen[i].wgt());
      }
    }
    //---------------------------------------------------------------

    // Do weight cancelation
    if (settings->regional_cancellation) {
      perform_regional_cancelation(next_gen);
    }

    // Calculate net positive and negative weight
    double W = 0.;
    double W_neg = 0.;
    double W_pos = 0.;
    int N_net = 0;
    int N_tot = 0;
    int N_pos = 0;
    int N_neg = 0;
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < next_gen.size(); i++) {
      if (next_gen[i].wgt() > 0.) {
#pragma omp atomic
        W_pos += next_gen[i].wgt();
#pragma omp atomic
        N_pos++;
      } else {
#pragma omp atomic
        W_neg -= next_gen[i].wgt();
#pragma omp atomic
        N_neg++;
      }
    }
    W = W_pos - W_neg;
    N_net = N_pos - N_neg;
    N_tot = N_pos + N_neg;

    // Re-Normalize particle weights
    double w_per_part = static_cast<double>(settings->nparticles) / W;
    W_neg = 0.;
    W_pos = 0.;
    W = 0.;
#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < next_gen.size(); i++) {
      next_gen[i].set_weight(next_gen[i].wgt() * w_per_part);

#pragma omp atomic
      W += next_gen[i].wgt();

      if (next_gen[i].wgt() > 0.) {
#pragma omp atomic
        W_pos += next_gen[i].wgt();
      } else {
#pragma omp atomic
        W_neg -= next_gen[i].wgt();
      }

      //---------------------------------------------------------------
      // Do all Post-Cancelation entropy calculations
      if (t_pre_entropy) {
        p_post_entropy->add_point(next_gen[i].r(), next_gen[i].wgt());
        n_post_entropy->add_point(next_gen[i].r(), next_gen[i].wgt());
        t_post_entropy->add_point(next_gen[i].r(), next_gen[i].wgt());
      }
      //---------------------------------------------------------------
    }

    double Wtt = W_pos + W_neg;
    int Wnet = std::round(W);
    int Wtot = std::round(Wtt);
    int Wpos = std::round(W_pos);
    int Wneg = std::round(W_neg);

    // Clear and switch particle banks
    bank.clear();
    bank = next_gen;
    next_gen.clear();

    // Output
    generation_output(g, N_net, N_tot, N_pos, N_neg, Wnet, Wtot, Wpos, Wneg);

    // Check if signal has been sent after generation keff has been
    // recorded, and cancellation has occured. Otherwize source file
    // will be larger than necessary, and wrong number of gens will be
    // in output file based on number averaged for tallies.
    if (signaled) premature_kill();

    // Zero all entropy bins
    if (t_pre_entropy) {
      p_pre_entropy->zero();
      n_pre_entropy->zero();
      t_pre_entropy->zero();

      p_post_entropy->zero();
      n_post_entropy->zero();
      t_post_entropy->zero();
    }

    // Once ignored generations are finished, mark as true to start
    // doing tallies
    if (g == settings->nignored) settings->converged = true;
  }

  Output::instance()->write("\n");

  // Write flux file
  if (settings->converged) {
    tallies->write_flux(settings->flux_file_name);
    tallies->write_power(settings->power_file_name);
  }

  // Check to write source
  if (settings->save_source) write_source(bank, settings->source_file_name);
}

void PowerIterator::generation_output(int gen, int Nn, int Nt, int Npos,
                                      int Nneg, int Wn, int Wt, int Wpos,
                                      int Wneg) {
  std::shared_ptr<Output> out = Output::instance();
  std::stringstream output;

  // Print header again every 300 generations
  if (((gen - 1) % 300) == 0 && (gen - 1 != settings->nignored)) print_header();

  if (gen <= settings->nignored + 1) {
    output << std::fixed << std::setw(5) << std::setfill(' ') << gen << "   "
           << std::setprecision(5) << tallies->keff();
    if (t_pre_entropy) {
      output << "                         ";
      output << std::setw(10) << std::scientific << std::setprecision(4)
             << p_pre_entropy->calculate_entropy() << "   ";
      output << n_pre_entropy->calculate_entropy() << "   ";
      output << t_pre_entropy->calculate_entropy() << "   ";

      output << p_post_entropy->calculate_entropy() << "   ";
      output << n_post_entropy->calculate_entropy() << "   ";
      output << t_post_entropy->calculate_entropy() << "   ";
    } else
      output << "                         ";
    output << std::fixed << std::setw(7) << std::setfill(' ') << Nn << "   ";
    output << std::fixed << std::setw(7) << std::setfill(' ') << Nt << "   ";
    output << std::fixed << std::setw(7) << std::setfill(' ') << Npos << "   ";
    output << std::fixed << std::setw(7) << std::setfill(' ') << Nneg << "   ";
    output << std::fixed << std::setw(7) << std::setfill(' ') << Wn << "   ";
    output << std::fixed << std::setw(7) << std::setfill(' ') << Wt << "   ";
    output << std::fixed << std::setw(7) << std::setfill(' ') << Wpos << "   ";
    output << std::fixed << std::setw(7) << std::setfill(' ') << Wneg
           << "   \n";
    out->write(output.str());
    output.str(std::string());
    if (gen == settings->nignored) {
      print_header();
      /*if(t_pre_entropy)
        out->write("
      ----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
      else
        out->write("
      ----------------------------------------------------------------------------\n");*/
    }
  } else {
    output << std::fixed << std::setw(5) << std::setfill(' ') << gen << "   ";
    output << std::setprecision(5) << tallies->keff() << "   ";
    output << tallies->kavg() << " +/- " << tallies->kerr() << "   ";
    if (t_pre_entropy) {
      output << std::setw(10) << std::scientific << std::setprecision(4)
             << p_pre_entropy->calculate_entropy() << "   ";
      output << n_pre_entropy->calculate_entropy() << "   ";
      output << t_pre_entropy->calculate_entropy() << "   ";

      output << p_post_entropy->calculate_entropy() << "   ";
      output << n_post_entropy->calculate_entropy() << "   ";
      output << t_post_entropy->calculate_entropy() << "   ";
    }
    output << std::fixed << std::setw(7) << std::setfill(' ') << Nn << "   ";
    output << std::fixed << std::setw(7) << std::setfill(' ') << Nt << "   ";
    output << std::fixed << std::setw(7) << std::setfill(' ') << Npos << "   ";
    output << std::fixed << std::setw(7) << std::setfill(' ') << Nneg << "   ";
    output << std::fixed << std::setw(7) << std::setfill(' ') << Wn << "   ";
    output << std::fixed << std::setw(7) << std::setfill(' ') << Wt << "   ";
    output << std::fixed << std::setw(7) << std::setfill(' ') << Wpos << "   ";
    output << std::fixed << std::setw(7) << std::setfill(' ') << Wneg
           << "   \n";
    out->write(output.str());
    output.str(std::string());
  }
}

void PowerIterator::write_source(std::vector<Particle>& bank,
                                 std::string source_fname) {
  Output::instance()->write(" Writing source distribution file...\n");

  // Make an NDArray to contain all particles info first
  NDArray<double> source({bank.size(), 5});

  // Add all particles to the array
  for (size_t i = 0; i < bank.size(); i++) {
    source[i * 5 + 0] = bank[i].r().x();
    source[i * 5 + 1] = bank[i].r().y();
    source[i * 5 + 2] = bank[i].r().z();
    source[i * 5 + 3] = static_cast<double>(bank[i].E());
    source[i * 5 + 4] = bank[i].wgt();
  }

  source.save(source_fname + ".npy");
}

void PowerIterator::premature_kill() {
  // See if the user really wants to kill the program
  std::shared_ptr<Output> out = Output::instance();
  std::string response;
  std::cout << "\n Do you really want to stop the simulation ? (y/N) => ";
  std::cin >> response;

  if (response == "y" || response == "Y") {
    out->write("\n Simulation has been stopped by user. Cleaning up...\n");
    // Looks like they really want to kill it.
    // Save everthing for them, and abort.

    // Only save flux is sim had already reached active generations.
    // Otherwise there will be no flux to write
    if (settings->converged) {
      tallies->write_flux(settings->flux_file_name);
      tallies->write_power(settings->power_file_name);
    }

    // Only save source if it was requested in the input file
    if (settings->save_source) write_source(bank, settings->source_file_name);

    // Finally, we may gracefully kill the program
    std::exit(0);
  }

  // They don't really want to kill it. Reset flag
  signaled = false;
  std::cout << "\n";
}

std::string PowerIterator::remove(char c, std::string line) {
  std::string out_line = "";
  for (size_t i = 0; i < line.size(); i++) {
    if (line[i] != c) out_line += line[i];
  }

  return out_line;
}

std::vector<std::string> PowerIterator::split(std::string line) {
  std::vector<std::string> output;
  std::string element = "";

  for (size_t i = 0; i < line.size(); i++) {
    if (line[i] != ',')
      element += line[i];
    else {
      if (element.size() != 0) {
        output.push_back(element);
        element.clear();
      }
    }
  }
  if (element.size() > 0) output.push_back(element);
  return output;
}

void PowerIterator::perform_regional_cancelation(
    std::vector<Particle>& next_gen) {
  size_t n_lost_boys = 0;

  double Xl = geometry::cancel_bins_low.x();
  double Yl = geometry::cancel_bins_low.y();
  double Zl = geometry::cancel_bins_low.z();

  double Px = geometry::cancel_bins_pitch[0];
  double Py = geometry::cancel_bins_pitch[1];
  double Pz = geometry::cancel_bins_pitch[2];

  size_t Nx = geometry::cancel_bins_shape[0];
  size_t Ny = geometry::cancel_bins_shape[1];
  size_t Nz = geometry::cancel_bins_shape[2];

  // CANNOT RUN IN PARALLEL DUE TO ADDING PARTICLES
  for (size_t i = 0; i < next_gen.size(); i++) {
    // Calculate indicies for particle
    int nx = static_cast<int>(std::floor((next_gen[i].r().x() - Xl) / Px));
    int ny = static_cast<int>(std::floor((next_gen[i].r().y() - Yl) / Py));
    int nz = static_cast<int>(std::floor((next_gen[i].r().z() - Zl) / Pz));

    // Only add particle if the bin is inside the mesh
    if (nx >= 0 && nx < static_cast<int>(Nx) && ny >= 0 &&
        ny < static_cast<int>(Ny) && nz >= 0 && nz < static_cast<int>(Nz)) {
      CancelBinKey binKey{static_cast<size_t>(nx), static_cast<size_t>(ny),
                          static_cast<size_t>(nz)};

      double Xl_bin =
          geometry::cancel_bins_low.x() + nx * geometry::cancel_bins_pitch[0];
      double Xh_bin = Xl_bin + geometry::cancel_bins_pitch[0];
      double Yl_bin =
          geometry::cancel_bins_low.y() + ny * geometry::cancel_bins_pitch[1];
      double Yh_bin = Yl_bin + geometry::cancel_bins_pitch[1];
      double Zl_bin =
          geometry::cancel_bins_low.z() + nz * geometry::cancel_bins_pitch[2];
      double Zh_bin = Zl_bin + geometry::cancel_bins_pitch[2];

      // Get bin center for material
      Position binCenter(0.5 * (Xl_bin + Xh_bin), 0.5 * (Yl_bin + Yh_bin),
                         0.5 * (Zl_bin + Zh_bin));

      // Get bin cell/material
      std::shared_ptr<Cell> binCell =
          geometry::get_cell(binCenter, {1., 0., 0.}, 0);
      std::shared_ptr<Material> binMat = binCell->material();

      // Get material at particle location
      std::shared_ptr<Cell> prtCell =
          geometry::get_cell(next_gen[i].r(), next_gen[i].u(), 0);
      std::shared_ptr<Material> prtMat = prtCell->material();

      // Make sure binMat and prtMat are the same material
      if (binMat->id() == prtMat->id()) {
        // Check if bin exists
        if (geometry::cancel_bins.find(binKey) != geometry::cancel_bins.end()) {
          // Bin exists, add particle
          geometry::cancel_bins[binKey].add_particle(
              next_gen[i], Xl_bin, Xh_bin, Yl_bin, Yh_bin, Zl_bin, Zh_bin);
        } else {
          // Create bin
          geometry::cancel_bins[binKey] = CancelBin(binMat);

          // Add particle
          geometry::cancel_bins[binKey].add_particle(
              next_gen[i], Xl_bin, Xh_bin, Yl_bin, Yh_bin, Zl_bin, Zh_bin);
        }
      } else
        n_lost_boys += 1;
    } else
      n_lost_boys += 1;

  }  // End distribute all particle to CancelBins

  if (n_lost_boys > 0)
    std::cout << " There are " << n_lost_boys
              << " particles with no cancellation bin.\n";

  // Do cancelation for each bin. Can run in parallel. First get keys
  std::vector<CancelBinKey> keys;
  keys.reserve(geometry::cancel_bins.size());
  for (auto it = geometry::cancel_bins.begin();
       it != geometry::cancel_bins.end(); it++) {
    keys.push_back(it->first);
  }

#pragma omp parallel
  {
#pragma omp for schedule(static)
    for (size_t k = 0; k < keys.size(); k++) {
      CancelBinKey key = keys[k];

      double Xl_bin = geometry::cancel_bins_low.x() +
                      key.i * geometry::cancel_bins_pitch[0];
      double Xh_bin = Xl_bin + geometry::cancel_bins_pitch[0];
      double Yl_bin = geometry::cancel_bins_low.y() +
                      key.j * geometry::cancel_bins_pitch[1];
      double Yh_bin = Yl_bin + geometry::cancel_bins_pitch[1];
      double Zl_bin = geometry::cancel_bins_low.z() +
                      key.k * geometry::cancel_bins_pitch[2];
      double Zh_bin = Zl_bin + geometry::cancel_bins_pitch[2];

      geometry::cancel_bins[key].perform_cancelation(Xl_bin, Xh_bin, Yl_bin,
                                                     Yh_bin, Zl_bin, Zh_bin);
    }
  }
  // All particles which were placed into a cancellation bin from next_gen
  // now have modified weights.

  // Get uniform particles CANNOT DO IN PARALLEL
  auto rng = rngs[0];
  for (size_t k = 0; k < keys.size(); k++) {
    CancelBinKey key = keys[k];

    double Xl_bin =
        geometry::cancel_bins_low.x() + key.i * geometry::cancel_bins_pitch[0];
    double Xh_bin = Xl_bin + geometry::cancel_bins_pitch[0];
    double Yl_bin =
        geometry::cancel_bins_low.y() + key.j * geometry::cancel_bins_pitch[1];
    double Yh_bin = Yl_bin + geometry::cancel_bins_pitch[1];
    double Zl_bin =
        geometry::cancel_bins_low.z() + key.k * geometry::cancel_bins_pitch[2];
    double Zh_bin = Zl_bin + geometry::cancel_bins_pitch[2];

    auto tmp = geometry::cancel_bins[key].get_uniform_particles(
        rng, Xl_bin, Xh_bin, Yl_bin, Yh_bin, Zl_bin, Zh_bin);
    next_gen.insert(next_gen.end(), tmp.begin(), tmp.end());
  }
}
