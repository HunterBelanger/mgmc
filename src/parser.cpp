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
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <geometry/cell_universe.hpp>
#include <geometry/geometry.hpp>
#include <geometry/hex_lattice.hpp>
#include <geometry/lattice_universe.hpp>
#include <geometry/rect_lattice.hpp>
#include <geometry/surfaces/all_surfaces.hpp>
#include <iomanip>
#include <ios>
#include <materials/nuclide.hpp>
#include <simulation/carter_tracker.hpp>
#include <simulation/delta_tracker.hpp>
#include <simulation/entropy.hpp>
#include <simulation/fixed_source.hpp>
#include <simulation/mesh_tally.hpp>
#include <simulation/noise.hpp>
#include <simulation/power_iterator.hpp>
#include <simulation/surface_tracker.hpp>
#include <sstream>
#include <string>
#include <utils/error.hpp>
#include <utils/output.hpp>
#include <utils/parser.hpp>
#include <utils/settings.hpp>
#include <utils/timer.hpp>

//===========================================================================
// Initialize Maps from id to index
std::map<uint32_t, size_t> surface_id_to_indx;
std::map<uint32_t, size_t> cell_id_to_indx;
std::map<uint32_t, size_t> universe_id_to_indx;
std::map<uint32_t, size_t> lattice_id_to_indx;

//===========================================================================
// Objects to build Simulation
std::vector<std::shared_ptr<Source>> sources;
NoiseMaker noise_maker;
std::shared_ptr<Tallies> tallies = nullptr;
std::shared_ptr<Transporter> transporter = nullptr;
std::shared_ptr<Simulation> simulation = nullptr;
std::shared_ptr<Cancelator> cancelator = nullptr;
std::string xspath;

void parse_input_file(std::string fname) {
  // Start parsing timer
  Timer parsing_timer;
  parsing_timer.start();

  Output::instance()->write(" Reading input file.\n");

  // Open input file
  YAML::Node input = YAML::LoadFile(fname);

  make_settings(input);

  Output::instance()->write(" Constructing simulation.\n");
  make_materials(input);

  make_geometry(input);

  make_tallies(input);

  make_transporter();

  if (settings::regional_cancellation ||
      settings::regional_cancellation_noise) {
    make_cancellation_bins(input);
  }

  make_sources(input);

  make_noise_sources(input);

  make_simulation();

  // Only parse entropy mesh if there is an entropy entry
  if (input["entropy"] && input["entropy"].IsMap())
    make_entropy_mesh(input["entropy"]);

  // End parsing timer
  parsing_timer.stop();
  Output::instance()->write(
      " Time to Parse Input : " + std::to_string(parsing_timer.elapsed_time()) +
      " seconds.\n");
}

void make_materials(YAML::Node input, bool plotting_mode) {
  // Make an initially Null xsdir node
  YAML::Node xsdir;

  // Parse materials
  if (input["materials"] && input["materials"].IsSequence()) {
    // Go through all materials
    for (size_t s = 0; s < input["materials"].size(); s++) {
      make_material(input["materials"][s], xsdir, plotting_mode);
    }

  } else {
    // If materials doesn't exist and isn't a sequence, kill program
    std::string mssg = "No materials are provided in input file.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Once all materials have been read in, we need to go find the max
  // and min energy grid values. To do this, we loop through all
  // nuclides in the problem.
  if (plotting_mode) return;

  for (const auto &mat : materials) {
    double emax = mat.second->max_energy();
    double emin = mat.second->min_energy();

    if (emin > settings::min_energy) {
      settings::min_energy = emin;
    }

    if (emax < settings::max_energy) {
      settings::max_energy = emax;
    }
  }

  std::stringstream mssg;
  mssg << " Min energy = " << std::setprecision(3) << std::scientific
       << settings::min_energy << " MeV.\n";
  mssg << " Max energy = " << settings::max_energy << " MeV.\n";

  Output::instance()->write(mssg.str());
}

void make_geometry(YAML::Node input) {
  std::shared_ptr<Output> out = Output::instance();

  // Parse Surfaces
  if (input["surfaces"] && input["surfaces"].IsSequence()) {
    // Go through all surfaces
    for (size_t s = 0; s < input["surfaces"].size(); s++) {
      make_surface(input["surfaces"][s]);
    }

  } else {
    // If surfaces doesn't exist and isn't a sequence, kill program
    std::string mssg = "No surfaces are provided in input file.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Parse Cells
  if (input["cells"] && input["cells"].IsSequence()) {
    // Go through all cells
    for (size_t c = 0; c < input["cells"].size(); c++) {
      make_cell(input["cells"][c]);
    }

  } else {
    // If cells doesn't exist and isn't a sequence, kill program
    std::string mssg = "No cells are provided in input file.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Parse all universes
  if (input["universes"] && input["universes"].IsSequence()) {
    // Go through all cells
    for (size_t c = 0; c < input["universes"].size(); c++) {
      make_universe(input["universes"][c], input);
    }

  } else {
    // If cells doesn't exist and isn't a sequence, kill program
    std::string mssg = "No universes are provided in input file.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Parse root universe
  if (input["root-universe"] && input["root-universe"].IsScalar()) {
    uint32_t root_id = input["root-universe"].as<uint32_t>();

    // Make sure it can be found
    if (universe_id_to_indx.find(root_id) == universe_id_to_indx.end()) {
      std::string mssg = "Root-Universe id was not found.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    geometry::root_universe = geometry::universes[universe_id_to_indx[root_id]];
  } else {
    // If doesn't exist and isn't a scalar, kill program
    std::string mssg = "No root-universe is provided in input file.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // If there is only one universe, MGMC allows for the user to
  // provide a neighbors list, which is a dictionary of entries, one for
  // each cell, and each entry is a sequence of the IDs of all it's
  // neighboring cells. This speeds up transport in certain cases.
  if (geometry::universes.size() == 1 && input["neighbors"] &&
      input["neighbors"].IsSequence()) {
    out->write(" Reading neighbors lists.\n");

    auto neighbors_array =
        input["neighbors"].as<std::vector<std::vector<uint32_t>>>();

    // Go through all cells in the geometry
    for (const auto &cell : geometry::cells) {
      // Get reference to vector of all neighbor IDs
      const std::vector<uint32_t> &neighbors = neighbors_array[cell->id() - 1];

      for (const auto &neighbor : neighbors) {
        // Check if the neighbor doesn't exist
        if (cell_id_to_indx.find(neighbor) == cell_id_to_indx.end()) {
          // We couldn't find a cell with that ID, so we skip trying to
          // add this neighbor.
          continue;
        }

        // Get the cell index for the cell ID
        auto nbr_indx = cell_id_to_indx[neighbor];

        // Add neighbor to cell
        cell->add_neighbor(geometry::cells[nbr_indx]);
      }
    }
  } else if (geometry::universes.size() == 1 && input["neighbors"]) {
    std::string mssg = "Invalid neighbors format.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // If there is a cell search mesh for the stochastic geometry, we now
  // read that too.
  if (geometry::universes.size() == 1 && input["cell-search-mesh"] &&
      input["cell-search-mesh"].IsMap()) {
    out->write(" Reading cell search mesh.\n");
    geometry::cell_search_mesh =
        make_cell_search_mesh(input["cell-search-mesh"]);
  } else if (geometry::universes.size() == 1 && input["cell-search-mesh"]) {
    std::string mssg = "Invalid cell-search-mesh format.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
}

void make_surface(YAML::Node surface_node) {
  // Try to get type
  std::string surf_type;
  if (surface_node["type"] && surface_node["type"].IsScalar())
    surf_type = surface_node["type"].as<std::string>();
  else {
    // Error, all surfaces must have a type
    std::string mssg = "Surface is missing \"type\" attribute.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Call appropriate function to build pointer to surface
  std::shared_ptr<Surface> surf_pntr = nullptr;
  if (surf_type == "xplane") {
    surf_pntr = make_xplane(surface_node);
  } else if (surf_type == "yplane") {
    surf_pntr = make_yplane(surface_node);
  } else if (surf_type == "zplane") {
    surf_pntr = make_zplane(surface_node);
  } else if (surf_type == "plane") {
    surf_pntr = make_plane(surface_node);
  } else if (surf_type == "xcylinder") {
    surf_pntr = make_xcylinder(surface_node);
  } else if (surf_type == "ycylinder") {
    surf_pntr = make_ycylinder(surface_node);
  } else if (surf_type == "zcylinder") {
    surf_pntr = make_zcylinder(surface_node);
  } else if (surf_type == "sphere") {
    surf_pntr = make_sphere(surface_node);
  } else {
    // Error, unknown surface type
    std::string mssg = "Surface type \"" + surf_type + "\" is unknown.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Add surface ID to map of surface indicies
  if (surface_id_to_indx.find(surf_pntr->id()) != surface_id_to_indx.end()) {
    // ID already exists
    std::string mssg = "The surface ID " + std::to_string(surf_pntr->id()) +
                       " appears multiple times.";
    fatal_error(mssg, __FILE__, __LINE__);
  } else {
    surface_id_to_indx[surf_pntr->id()] = geometry::surfaces.size();
    geometry::surfaces.push_back(surf_pntr);
  }

  // Add index to boundaries vector, which contains indicies of surfaces
  // which have critical boundary conditions
  if (surf_pntr->boundary() != BoundaryType::Normal) {
    geometry::boundaries.push_back(surface_id_to_indx[surf_pntr->id()]);
  }
}

void make_universe(YAML::Node uni_node, YAML::Node input) {
  uint32_t id;
  if (uni_node["id"] && uni_node["id"].IsScalar()) {
    id = uni_node["id"].as<uint32_t>();

    // Make sure cell can't be found
    if (universe_id_to_indx.find(id) == universe_id_to_indx.end()) {
      if (uni_node["cells"] && uni_node["cells"].IsSequence()) {
        make_cell_universe(uni_node);
      } else if (uni_node["lattice"] && uni_node["lattice"].IsScalar()) {
        make_lattice_universe(uni_node, input);
      } else {
        std::string mssg = "Invalid universe definition.";
        fatal_error(mssg, __FILE__, __LINE__);
      }
    }

  } else {
    std::string mssg = "Universe must have a valid id.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
}

void find_universe(YAML::Node input, uint32_t id) {
  // Iterate through universes
  if (universe_id_to_indx.find(id) == universe_id_to_indx.end()) {
    bool found = false;
    for (size_t u = 0; u < input["universes"].size(); u++) {
      if (input["universes"][u]["id"] &&
          input["universes"][u]["id"].IsScalar()) {
        uint32_t u_id = input["universes"][u]["id"].as<uint32_t>();
        if (u_id == id) {
          make_universe(input["universes"][u], input);
          found = true;
          break;
        }
      } else {
        std::string mssg = "Universes must have a valid id.";
        fatal_error(mssg, __FILE__, __LINE__);
      }
    }

    if (found == false) {
      std::string mssg =
          "Could not find universe with id " + std::to_string(id) + ".";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  }
}

void make_settings(YAML::Node input) {
  if (input["settings"] && input["settings"].IsMap()) {
    const auto &settnode = input["settings"];

    // Get simulation type
    std::string sim_type;
    if (settnode["simulation"] && settnode["simulation"].IsScalar()) {
      sim_type = settnode["simulation"].as<std::string>();
    } else {
      std::string mssg = "No simulation type provided.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    if (sim_type == "k-eigenvalue") {
      settings::mode = settings::SimulationMode::K_EIGENVALUE;
    } else if (sim_type == "fixed-source") {
      settings::mode = settings::SimulationMode::FIXED_SOURCE;
    } else if (sim_type == "noise") {
      settings::mode = settings::SimulationMode::NOISE;
    } else {
      std::string mssg = "Unknown simulation type " + sim_type + ".";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Get transport method
    if (settnode["transport"] && settnode["transport"].IsScalar()) {
      if (settnode["transport"].as<std::string>() == "delta-tracking") {
        settings::tracking = settings::TrackingMode::DELTA_TRACKING;
      } else if (settnode["transport"].as<std::string>() ==
                 "surface-tracking") {
        settings::tracking = settings::TrackingMode::SURFACE_TRACKING;
      } else if (settnode["transport"].as<std::string>() == "carter-tracking") {
        settings::tracking = settings::TrackingMode::CARTER_TRACKING;

        // Read the sampling-xs entry in the input file
        if (!input["sampling-xs"] || !input["sampling-xs"].IsSequence()) {
          std::string mssg =
              "Must provide the \"sampling-xs\" vector to use Carter Tracking.";
          fatal_error(mssg, __FILE__, __LINE__);
        }

        settings::sample_xs_ratio =
            input["sampling-xs"].as<std::vector<double>>();

        // Make sure all ratios are positive
        for (const auto &v : settings::sample_xs_ratio) {
          if (v <= 0.) {
            std::string mssg = "Sampling XS ratios must be > 0.";
            fatal_error(mssg, __FILE__, __LINE__);
          }
        }
      } else {
        std::string mssg = "Invalid tracking method " +
                           settnode["transport"].as<std::string>() + ".";
        fatal_error(mssg, __FILE__, __LINE__);
      }
    } else {
      // By default, use surface tracking
      settings::tracking = settings::TrackingMode::SURFACE_TRACKING;
    }

    // Get energy mode type
    if (settnode["energy-mode"] && settnode["energy-mode"].IsScalar()) {
      if (settnode["energy-mode"].as<std::string>() == "multi-group") {
        settings::energy_mode = settings::EnergyMode::MG;
        Output::instance()->write(" Running in Multi-Group mode.\n");
      } else if (settnode["energy-mode"].as<std::string>() ==
                 "continuous-energy") {
        std::string mssg = "Continuous-Energy Mode is not supported.";
        fatal_error(mssg, __FILE__, __LINE__);
      } else {
        std::string mssg = "Invalid energy mode ";
        mssg += settnode["energy-mode"].as<std::string>() + ".";
        fatal_error(mssg, __FILE__, __LINE__);
      }
    } else {
      // MG by default
      settings::energy_mode = settings::EnergyMode::MG;
      Output::instance()->write(" Running in Multi-Group mode.\n");
    }

    // If we are multi-group, get number of groups
    if (settings::energy_mode == settings::EnergyMode::MG) {
      if (!settnode["ngroups"] || !settnode["ngroups"].IsScalar()) {
        std::string mssg =
            "Number of groups for mutli-group mode not provided.";
        fatal_error(mssg, __FILE__, __LINE__);
      }
      int ngrps = settnode["ngroups"].as<int>();

      if (ngrps <= 0) {
        std::string mssg = "Number of groups may not be negative.";
        fatal_error(mssg, __FILE__, __LINE__);
      }

      settings::ngroups = static_cast<uint32_t>(ngrps);

      // Must also get the energy-bounds for the groups
      if (!settnode["energy-bounds"] ||
          !settnode["energy-bounds"].IsSequence()) {
        std::string mssg =
            "No energy-bounds entry found in settings for multi-group mode.";
        fatal_error(mssg, __FILE__, __LINE__);
      }

      if (settnode["energy-bounds"].size() != settings::ngroups + 1) {
        std::string mssg =
            "The number of energy-bounds must be equal to ngroups + 1.";
        fatal_error(mssg, __FILE__, __LINE__);
      }
      settings::energy_bounds =
          settnode["energy-bounds"].as<std::vector<double>>();

      // Check bounds
      if (!std::is_sorted(settings::energy_bounds.begin(),
                          settings::energy_bounds.end())) {
        std::string mssg =
            "The energy-bounds for multi-group mode are not sorted.";
        fatal_error(mssg, __FILE__, __LINE__);
      }

      // If we are using carter tracking, make sure we have the right number
      // of sampling xs ratios !
      if (settings::tracking == settings::TrackingMode::CARTER_TRACKING) {
        if (settings::sample_xs_ratio.size() != settings::ngroups) {
          std::string mssg =
              "The number of energy groups does not match the size of "
              "\"sampling-xs\".";
          fatal_error(mssg, __FILE__, __LINE__);
        }
      }
    }

    // Get number of particles
    if (settnode["nparticles"] && settnode["nparticles"].IsScalar()) {
      settings::nparticles = settnode["nparticles"].as<int>();
    } else {
      std::string mssg = "Number of particles not specified in settings.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Get number of generations
    if (settnode["ngenerations"] && settnode["ngenerations"].IsScalar()) {
      settings::ngenerations = settnode["ngenerations"].as<int>();
    } else {
      std::string mssg = "Number of generations not specified in settings.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Get number of ignored
    if (settnode["nignored"] && settnode["nignored"].IsScalar()) {
      settings::nignored = settnode["nignored"].as<int>();
    } else if (settnode["simulation"] &&
               settnode["simulation"].as<std::string>() == "k-eigenvalue") {
      std::string mssg =
          "Number of ignored generations not specified in settings.";
      fatal_error(mssg, __FILE__, __LINE__);
    } else {
      settings::nignored = 0;
    }

    // Get number of skips between noise batches.
    // Default is 10
    if (settnode["nskip"] && settnode["nskip"].IsScalar()) {
      settings::nskip = settnode["nskip"].as<int>();
    }

    // Get max run time
    if (settnode["max-run-time"] && settnode["max-run-time"].IsScalar()) {
      // Get runtime in minutes
      double max_time_in_mins = settnode["max-run-time"].as<double>();
      Output::instance()->write(
          " Max Run Time: " + std::to_string(max_time_in_mins) + " mins.\n");

      // Change minutes to seconds
      settings::max_time = max_time_in_mins * 60.;
    } else if (settnode["max-run-time"]) {
      std::string mssg = "Invalid max-run-time entry in settings.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Get the frequency and keff for noise simulations
    if (settings::mode == settings::SimulationMode::NOISE) {
      // Frequency
      if (!settnode["noise-angular-frequency"] ||
          !settnode["noise-angular-frequency"].IsScalar()) {
        std::string mssg =
            "No valid noise-angular-frequency in settings for noise "
            "simulation.";
        fatal_error(mssg, __FILE__, __LINE__);
      }
      settings::w_noise = settnode["noise-angular-frequency"].as<double>();

      Output::instance()->write(
          " Noise angular-frequency: " + std::to_string(settings::w_noise) +
          " radians / s.\n");

      // Keff
      if (!settnode["keff"] || !settnode["keff"].IsScalar()) {
        std::string mssg = "No valid keff in settings for noise simulation.";
        fatal_error(mssg, __FILE__, __LINE__);
      }
      settings::keff = settnode["keff"].as<double>();

      if (settings::keff <= 0.) {
        std::string mssg = "Noise keff must be greater than zero.";
        fatal_error(mssg, __FILE__, __LINE__);
      }

      Output::instance()->write(
          " Noise keff: " + std::to_string(settings::keff) + "\n");

      // Get the inner_generations option
      if (settnode["inner-generations"] &&
          settnode["inner-generations"].IsScalar()) {
        settings::inner_generations = settnode["inner-generations"].as<bool>();
      } else if (settnode["inner-generations"]) {
        std::string mssg =
            "Invalid \"inner-generations\" entry provided in settings.";
        fatal_error(mssg, __FILE__, __LINE__);
      }

      // Get the normalize_noise_source option
      if (settnode["normalize-noise-source"] &&
          settnode["normalize-noise-source"].IsScalar()) {
        settings::normalize_noise_source =
            settnode["normalize-noise-source"].as<bool>();
        if (settings::normalize_noise_source) {
          Output::instance()->write(" Normalizing noise source\n");
        }
      } else if (settnode["normalize-noise-source"]) {
        std::string mssg =
            "Invalid \"normalize-noise-source\" entry provided in settings.";
        fatal_error(mssg, __FILE__, __LINE__);
      }
    }

    // See if the user wants to see RNG stride warnings
    if (settnode["rng-stride-warnings"] &&
        settnode["rng-stride-warnings"].IsScalar()) {
      settings::rng_stride_warnings =
          settnode["rng-stride-warnings"].as<bool>();
    } else if (settnode["rng-stride-warnings"]) {
      std::string mssg =
          "Invalid \"rng-stride-warnings\" entry provided in settings.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Get name of source out file
    if (settnode["sourcefile"] && settnode["sourcefile"].IsScalar()) {
      settings::source_file_name = settnode["sourcefile"].as<std::string>();
      settings::save_source = true;
    }

    // Get name of source in file
    if (settnode["insource"] && settnode["insource"].IsScalar()) {
      settings::in_source_file_name = settnode["insource"].as<std::string>();
      settings::load_source_file = true;

      // Make sure source file exists
      std::ifstream source_fl(settings::in_source_file_name);
      if (!source_fl.good()) {
        std::string mssg = "Could not find source input file with name \"";
        mssg += settings::in_source_file_name + "\".";
        fatal_error(mssg, __FILE__, __LINE__);
      }
    }

    // Get seed for rng
    if (settnode["seed"] && settnode["seed"].IsScalar()) {
      settings::rng_seed = settnode["seed"].as<uint64_t>();
      Output::instance()->write(
          " RNG seed = " + std::to_string(settings::rng_seed) + " ...\n");
    }

    // Get stride for rng
    if (settnode["stride"] && settnode["stride"].IsScalar()) {
      settings::rng_stride = settnode["stride"].as<uint64_t>();
      Output::instance()->write(
          " RNG stride = " + std::to_string(settings::rng_seed) + " ...\n");
    }

    // Get cancellation settings
    settings::regional_cancellation = false;
    if (settnode["cancellation"] && settnode["cancellation"].IsScalar()) {
      if (settnode["cancellation"].as<bool>() == true) {
        settings::regional_cancellation = true;
      }
    }
    if (settings::regional_cancellation) {
      Output::instance()->write(" Using regional cancellation.\n");
    }

    // Get number of noise generations to cancel
    if (settnode["cancel-noise-gens"] &&
        settnode["cancel-noise-gens"].IsScalar()) {
      settings::n_cancel_noise_gens = settnode["cancel-noise-gens"].as<int>();
    }

    settings::regional_cancellation_noise = false;
    if (settnode["noise-cancellation"] &&
        settnode["noise-cancellation"].IsScalar()) {
      if (settnode["noise-cancellation"].as<bool>() == true) {
        settings::regional_cancellation_noise = true;
        Output::instance()->write(
            " Cancel noise gens : " +
            std::to_string(settings::n_cancel_noise_gens) + " generations.\n");
      }
    }

  } else {
    std::string mssg = "Not settings specified in input file.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
}

void make_tallies(YAML::Node input) {
  // Make base tallies object which is required
  tallies =
      std::make_shared<Tallies>(static_cast<double>(settings::nparticles));

  if (!input["tallies"]) {
    return;
  } else if (!input["tallies"].IsSequence()) {
    std::string mssg = "Tallies entry must be provided as a sequence.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Add all spatial mesh tallies to the tallies instance
  for (size_t t = 0; t < input["tallies"].size(); t++) {
    add_mesh_tally(*tallies, input["tallies"][t]);
  }

  if (settings::mode == settings::SimulationMode::NOISE) {
    tallies->set_keff(settings::keff);
  }
}

void make_transporter() {
  switch (settings::tracking) {
    case settings::TrackingMode::SURFACE_TRACKING:
      transporter = std::make_shared<SurfaceTracker>(tallies);
      Output::instance()->write(" Using Surface-Tracking.\n");
      break;

    case settings::TrackingMode::DELTA_TRACKING:
      transporter = std::make_shared<DeltaTracker>(tallies);
      Output::instance()->write(" Using Delta-Tracking.\n");
      break;

    case settings::TrackingMode::CARTER_TRACKING:
      transporter = std::make_shared<CarterTracker>(tallies);
      Output::instance()->write(" Using Carter-Tracking.\n");
      break;
  }
}

void make_cancellation_bins(YAML::Node input) {
  if (!input["cancelator"] || !input["cancelator"].IsMap()) {
    std::string mssg =
        "Regional cancelation is activated, but no cancelator entry is "
        "provided.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make sure we are using cancelation with a valid transport method !
  if (settings::mode == settings::SimulationMode::FIXED_SOURCE) {
    std::string mssg =
        "Cancellation may only be used with k-eigenvalue "
        "or noise problems.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make cancelator will check the cancelator type agains the tracking type.
  // This is because exact cancelators may only be used with DT, or Carter
  // Tracking, while the approximate cancelator can be used with any tracking
  // method.
  cancelator = make_cancelator(input["cancelator"]);
}

void make_sources(YAML::Node input) {
  if (input["sources"] && input["sources"].IsSequence()) {
    // Do all sources
    for (size_t s = 0; s < input["sources"].size(); s++) {
      sources.push_back(make_source(input["sources"][s]));
    }
  } else if (settings::mode == settings::SimulationMode::FIXED_SOURCE ||
             !settings::load_source_file) {
    // No source file given
    std::string mssg = "No source specified for problem.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
}

void make_noise_sources(YAML::Node input) {
  if (input["noise-sources"] && input["noise-sources"].IsSequence()) {
    // Do all sources
    for (size_t s = 0; s < input["noise-sources"].size(); s++) {
      noise_maker.add_noise_source(input["noise-sources"][s]);
    }
  } else if (settings::mode == settings::SimulationMode::NOISE) {
    // No source file given
    std::string mssg = "No valid noise-source entry for noise problem.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  if (noise_maker.num_noise_sources() == 0 &&
      settings::mode == settings::SimulationMode::NOISE) {
    std::string mssg = "No noise source specified for noise problem.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
}

void make_simulation() {
  switch (settings::mode) {
    case settings::SimulationMode::K_EIGENVALUE:
      if (!settings::regional_cancellation) {
        simulation =
            std::make_shared<PowerIterator>(tallies, transporter, sources);
      } else {
        simulation = std::make_shared<PowerIterator>(tallies, transporter,
                                                     sources, cancelator);
      }
      break;

    case settings::SimulationMode::FIXED_SOURCE:
      simulation = std::make_shared<FixedSource>(tallies, transporter, sources);
      break;

    case settings::SimulationMode::NOISE:
      if (!settings::regional_cancellation &&
          !settings::regional_cancellation_noise) {
        simulation =
            std::make_shared<Noise>(tallies, transporter, sources, noise_maker);
      } else {
        simulation = std::make_shared<Noise>(tallies, transporter, sources,
                                             cancelator, noise_maker);
      }
      break;
  }
}

void make_entropy_mesh(YAML::Node entropy) {
  // Get lower corner
  std::vector<double> low;
  if (entropy["low"] && entropy["low"].IsSequence() &&
      entropy["low"].size() == 3) {
    low = entropy["low"].as<std::vector<double>>();
  } else {
    fatal_error("No valid lower corner provided for entropy mesh.", __FILE__,
                __LINE__);
  }

  // Get upper corner
  std::vector<double> hi;
  if (entropy["hi"] && entropy["hi"].IsSequence() &&
      entropy["hi"].size() == 3) {
    hi = entropy["hi"].as<std::vector<double>>();
  } else {
    fatal_error("No valid upper corner provided for entropy mesh.", __FILE__,
                __LINE__);
  }

  // Get shape
  std::vector<uint32_t> shape;
  if (entropy["shape"] && entropy["shape"].IsSequence() &&
      entropy["shape"].size() == 3) {
    shape = entropy["shape"].as<std::vector<uint32_t>>();
  } else {
    fatal_error("No valid shape provided for entropy mesh.", __FILE__,
                __LINE__);
  }

  // Add entropy to simulation
  Position low_r(low[0], low[1], low[2]);
  Position hi_r(hi[0], hi[1], hi[2]);
  std::array<uint32_t, 3> shp{shape[0], shape[1], shape[2]};

  simulation->set_p_pre_entropy(
      std::make_shared<Entropy>(low_r, hi_r, shp, Entropy::Sign::Positive));
  simulation->set_n_pre_entropy(
      std::make_shared<Entropy>(low_r, hi_r, shp, Entropy::Sign::Negative));
  simulation->set_t_pre_entropy(
      std::make_shared<Entropy>(low_r, hi_r, shp, Entropy::Sign::Total));

  simulation->set_p_post_entropy(
      std::make_shared<Entropy>(low_r, hi_r, shp, Entropy::Sign::Positive));
  simulation->set_n_post_entropy(
      std::make_shared<Entropy>(low_r, hi_r, shp, Entropy::Sign::Negative));
  simulation->set_t_post_entropy(
      std::make_shared<Entropy>(low_r, hi_r, shp, Entropy::Sign::Total));
}
