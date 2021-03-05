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
#include <geometry/cell_universe.hpp>
#include <geometry/geometry.hpp>
#include <geometry/hex_lattice.hpp>
#include <geometry/lattice_universe.hpp>
#include <geometry/rect_lattice.hpp>
#include <geometry/surfaces/all_surfaces.hpp>
#include <materials/const_material.hpp>
#include <simulation/cancel_bin.hpp>
#include <simulation/carter_tracker.hpp>
#include <simulation/delta_tracker.hpp>
#include <simulation/entropy.hpp>
#include <simulation/power_iterator.hpp>
#include <simulation/surface_tracker.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>
#include <utils/parser.hpp>
#include <utils/timer.hpp>

//===========================================================================
// Initialize Maps from id to index
std::map<uint32_t, size_t> surface_id_to_indx;
std::map<uint32_t, size_t> cell_id_to_indx;
std::map<uint32_t, size_t> universe_id_to_indx;
std::map<uint32_t, size_t> lattice_id_to_indx;

//===========================================================================
// Objects to build Simulation
std::shared_ptr<Settings> settings = nullptr;
std::vector<std::shared_ptr<Source>> sources;
std::shared_ptr<Tallies> tallies = nullptr;
std::shared_ptr<Transporter> transporter = nullptr;
std::shared_ptr<Simulation> simulation = nullptr;
bool using_carter_tracking = false;

void parse_input_file(std::string fname) {
  // Start parsing timer
  Timer::instance()->start_parsing_timer();

  Output::instance()->write(" Reading input file...\n");

  // Open input file
  YAML::Node input = YAML::LoadFile(fname);

  // Read materials
  make_materials(input);

  // Read and build geometry
  make_geometry(input);

  make_settings(input);

  make_tallies(input);

  make_transporter(input);

  make_cancellation_bins(input);

  make_sources(input);

  make_simulation(input);

  // Only parse entropy mesh if there is an entropy entry
  if (input["entropy"] && input["entropy"].IsMap())
    make_entropy_mesh(input["entropy"]);

  // End parsing timer
  Timer::instance()->end_parsing_timer();
}

void make_materials(YAML::Node input) {
  // Parse Surfaces
  if (input["materials"] && input["materials"].IsSequence()) {
    // Go through all materials
    for (size_t s = 0; s < input["materials"].size(); s++) {
      make_const_material(input["materials"][s]);
    }

  } else {
    // If materials doesn't exist and isn't a sequence, kill program
    std::string mssg = "No materials are provided in input file.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
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
  // Initialize settings pointer, all members are public
  settings = std::make_shared<Settings>();

  if (input["settings"] && input["settings"].IsMap()) {
    // Get number of particles
    if (input["settings"]["nparticles"] &&
        input["settings"]["nparticles"].IsScalar()) {
      settings->nparticles = input["settings"]["nparticles"].as<int>();
    } else {
      std::string mssg = "Number of particles not specified in settings.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Get number of generations
    if (input["settings"]["ngenerations"] &&
        input["settings"]["ngenerations"].IsScalar()) {
      settings->ngenerations = input["settings"]["ngenerations"].as<int>();
    } else {
      std::string mssg = "Number of generations not specified in settings.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Get number of ignored
    if (input["settings"]["nignored"] &&
        input["settings"]["nignored"].IsScalar()) {
      settings->nignored = input["settings"]["nignored"].as<int>();
    } else {
      std::string mssg =
          "Number of ignored generations not specified in settings.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Get number of groups
    if (input["settings"]["ngroups"] &&
        input["settings"]["ngroups"].IsScalar()) {
      settings->ngroups = input["settings"]["ngroups"].as<int>();
    } else {
      std::string mssg = "Number of groups not specified in settings.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Get flux out file
    if (input["settings"]["fluxfile"] &&
        input["settings"]["fluxfile"].IsScalar()) {
      settings->flux_file_name =
          input["settings"]["fluxfile"].as<std::string>();
    } else {
      settings->flux_file_name = "flux";
    }

    // Get flux out file
    if (input["settings"]["powerfile"] &&
        input["settings"]["powerfile"].IsScalar()) {
      settings->power_file_name =
          input["settings"]["powerfile"].as<std::string>();
    } else {
      settings->power_file_name = "power";
    }

    // Get name of source out file
    if (input["settings"]["sourcefile"] &&
        input["settings"]["sourcefile"].IsScalar()) {
      settings->source_file_name =
          input["settings"]["sourcefile"].as<std::string>();
      settings->save_source = true;
    }

    // Get name of source in file
    if (input["settings"]["insource"] &&
        input["settings"]["insource"].IsScalar()) {
      settings->in_source_file_name =
          input["settings"]["insource"].as<std::string>();
      settings->load_source_file = true;

      // Make sure source file exists
      std::ifstream source_fl(settings->in_source_file_name);
      if (!source_fl.good()) {
        std::string mssg = "Could not find source input file with name \"";
        mssg += settings->in_source_file_name + "\".";
        fatal_error(mssg, __FILE__, __LINE__);
      }
    }

    // Get seed for rng
    if (input["settings"]["seed"] && input["settings"]["seed"].IsScalar()) {
      settings->rng_seed = input["settings"]["seed"].as<uint64_t>();
      Output::instance()->write(
          " RNG seed = " + std::to_string(settings->rng_seed) + " ...\n");
    }

    // Determine which RNG to use
    if (input["settings"]["rng"] && input["settings"]["rng"].IsScalar()) {
      std::string rng_type = input["settings"]["rng"].as<std::string>();

      if (rng_type == "MT" || rng_type == "mt") {
        Output::instance()->write(" Using Mersenne Twister as PRNG...\n");
        settings->use_pcg = false;
      } else if (rng_type != "PCG" && rng_type != "pcg") {
        std::string mssg =
            "Unknown random number generator identifier " + rng_type + ".";
        fatal_error(mssg, __FILE__, __LINE__);
      } else {
        Output::instance()->write(" Using PCG as PRNG...\n");
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
      std::make_shared<Tallies>(static_cast<double>(settings->nparticles));

  // Add flux tally if defined
  if (input["flux"] && input["flux"].IsMap()) {
    // Get shape
    std::vector<int> shape;
    if (input["flux"]["shape"] && input["flux"]["shape"].IsSequence() &&
        input["flux"]["shape"].size() == 3) {
      shape = input["flux"]["shape"].as<std::vector<int>>();
    } else {
      std::string mssg = "No shape provided for flux.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Get lower position
    std::vector<double> low_coords;
    if (input["flux"]["low"] && input["flux"]["low"].IsSequence() &&
        input["flux"]["low"].size() == 3) {
      low_coords = input["flux"]["low"].as<std::vector<double>>();
    } else {
      std::string mssg = "No lower corner provided for flux.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Get upper coordinates
    std::vector<double> hi_coords;
    if (input["flux"]["hi"] && input["flux"]["hi"].IsSequence() &&
        input["flux"]["hi"].size() == 3) {
      hi_coords = input["flux"]["hi"].as<std::vector<double>>();
    } else {
      std::string mssg = "No upper corner provided for flux.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    Position low(low_coords[0], low_coords[1], low_coords[2]);
    Position hi(hi_coords[0], hi_coords[1], hi_coords[2]);
    int nx = shape[0];
    int ny = shape[1];
    int nz = shape[2];
    int ng = settings->ngroups;

    tallies->set_flux_tally(
        std::make_unique<FluxTally>(low, hi, nx, ny, nz, ng));
  }

  // Add power tally if defined
  if (input["power"] && input["power"].IsMap()) {
    // Get shape
    std::vector<int> shape;
    if (input["power"]["shape"] && input["power"]["shape"].IsSequence() &&
        input["power"]["shape"].size() == 3) {
      shape = input["power"]["shape"].as<std::vector<int>>();
    } else {
      std::string mssg = "No shape provided for power.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Get lower position
    std::vector<double> low_coords;
    if (input["power"]["low"] && input["power"]["low"].IsSequence() &&
        input["power"]["low"].size() == 3) {
      low_coords = input["power"]["low"].as<std::vector<double>>();
    } else {
      std::string mssg = "No lower corner provided for power.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Get upper coordinates
    std::vector<double> hi_coords;
    if (input["power"]["hi"] && input["power"]["hi"].IsSequence() &&
        input["power"]["hi"].size() == 3) {
      hi_coords = input["power"]["hi"].as<std::vector<double>>();
    } else {
      std::string mssg = "No upper corner provided for power.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    Position low(low_coords[0], low_coords[1], low_coords[2]);
    Position hi(hi_coords[0], hi_coords[1], hi_coords[2]);
    int nx = shape[0];
    int ny = shape[1];
    int nz = shape[2];

    tallies->set_power_tally(std::make_unique<PowerTally>(low, hi, nx, ny, nz));
  }
}

void make_transporter(YAML::Node input) {
  if (input["settings"]["transport"] &&
      input["settings"]["transport"].IsScalar()) {
    if (input["settings"]["transport"].as<std::string>() == "delta-tracking") {
      transporter =
          std::make_shared<DeltaTracker>(tallies, settings, majorant_xs);
    } else if (input["settings"]["transport"].as<std::string>() ==
               "carter-tracking") {
      // Get sampling XS ratios for carter tracking
      std::vector<double> sampling_xs;
      if (input["sampling-xs"] && input["sampling-xs"].IsSequence() &&
          static_cast<int>(input["sampling-xs"].size()) == settings->ngroups) {
        sampling_xs = input["sampling-xs"].as<std::vector<double>>();
      } else {
        std::string mssg = "No sampling cross section ratios provided.";
        fatal_error(mssg, __FILE__, __LINE__);
      }

      // Multiply all ratios by true majorant, already calculated
      for (size_t i = 0; i < majorant_xs.size(); i++)
        sampling_xs[i] *= majorant_xs[i];

      transporter =
          std::make_shared<CarterTracker>(tallies, settings, sampling_xs);
      using_carter_tracking = true;
    } else if (input["settings"]["transport"].as<std::string>() ==
               "surface-tracking") {
      transporter = std::make_shared<SurfaceTracker>(tallies, settings);
    } else {
      std::string mssg = "Invalid tracking method " +
                         input["settings"]["transport"].as<std::string>() + ".";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  } else {
    // Use delta-tracking as default
    transporter =
        std::make_shared<DeltaTracker>(tallies, settings, majorant_xs);
  }
}

void make_cancellation_bins(YAML::Node input) {
  // First check seetings and see if using regional cancelation
  if (input["settings"]["cancellation"] &&
      input["settings"]["cancellation"].IsScalar()) {
    if (input["settings"]["cancellation"].as<bool>() == true) {
      settings->regional_cancellation = true;

      // Get regional cancelation bins info

      // Get shape
      std::vector<size_t> shape;
      if (input["cancellation-bins"]["shape"] &&
          input["cancellation-bins"]["shape"].IsSequence() &&
          input["cancellation-bins"]["shape"].size() == 3) {
        shape = input["cancellation-bins"]["shape"].as<std::vector<size_t>>();
      } else {
        std::string mssg = "No shape given for cancellation bins.";
        fatal_error(mssg, __FILE__, __LINE__);
      }

      // Get low corner
      std::vector<double> low;
      if (input["cancellation-bins"]["low"] &&
          input["cancellation-bins"]["low"].IsSequence() &&
          input["cancellation-bins"]["low"].size() == 3) {
        low = input["cancellation-bins"]["low"].as<std::vector<double>>();
      } else {
        std::string mssg = "No low given for cancellation bins.";
        fatal_error(mssg, __FILE__, __LINE__);
      }

      // Get hi corner
      std::vector<double> hi;
      if (input["cancellation-bins"]["hi"] &&
          input["cancellation-bins"]["hi"].IsSequence() &&
          input["cancellation-bins"]["hi"].size() == 3) {
        hi = input["cancellation-bins"]["hi"].as<std::vector<double>>();
      } else {
        std::string mssg = "No hi given for cancellation bins.";
        fatal_error(mssg, __FILE__, __LINE__);
      }

      // Calculate pitch in z direction based on hi and low corners and shape
      double px = (hi[0] - low[0]) / static_cast<double>(shape[0]);
      double py = (hi[1] - low[1]) / static_cast<double>(shape[1]);
      double pz = (hi[2] - low[2]) / static_cast<double>(shape[2]);

      // Set the global parameters in geometry namespace
      geometry::cancel_bins_low = {low[0], low[1], low[2]};
      geometry::cancel_bins_shape = {shape[0], shape[1], shape[2]};
      geometry::cancel_bins_pitch = {px, py, pz};
    }
  }
}

void make_sources(YAML::Node input) {
  if (input["sources"] && input["sources"].IsSequence()) {
    // Do all sources
    for (size_t s = 0; s < input["sources"].size(); s++) {
      // Get weight
      double weight = 1.;
      if (input["sources"][s]["weight"] &&
          input["sources"][s]["weight"].IsScalar()) {
        weight = input["sources"][s]["weight"].as<double>();
      } else {
        std::string mssg = "No weight given to source.";
        fatal_error(mssg, __FILE__, __LINE__);
      }

      // Get chi
      std::vector<double> chi;
      if (input["sources"][s]["chi"] &&
          input["sources"][s]["chi"].IsSequence() &&
          static_cast<int>(input["sources"][s]["chi"].size()) <=
              settings->ngroups) {
        chi = input["sources"][s]["chi"].as<std::vector<double>>();
      } else {
        std::string mssg = "No chi given to source.";
        fatal_error(mssg, __FILE__, __LINE__);
      }

      // Get type
      std::string type;
      if (input["sources"][s]["type"] &&
          input["sources"][s]["type"].IsScalar()) {
        type = input["sources"][s]["type"].as<std::string>();
      } else {
        std::string mssg = "No type given to source.";
        fatal_error(mssg, __FILE__, __LINE__);
      }

      // Construct source and add to vector
      if (type == "box") {
        // Get low corner
        std::vector<double> low_coords;
        if (input["sources"][s]["low"] &&
            input["sources"][s]["low"].IsSequence() &&
            input["sources"][s]["low"].size() == 3) {
          low_coords = input["sources"][s]["low"].as<std::vector<double>>();
        } else {
          std::string mssg = "Now lower coordinates given for box source.";
          fatal_error(mssg, __FILE__, __LINE__);
        }

        // Get hi corner
        std::vector<double> hi_coords;
        if (input["sources"][s]["hi"] &&
            input["sources"][s]["hi"].IsSequence() &&
            input["sources"][s]["hi"].size() == 3) {
          hi_coords = input["sources"][s]["hi"].as<std::vector<double>>();
        } else {
          std::string mssg = "Now lower coordinates given for box source.";
          fatal_error(mssg, __FILE__, __LINE__);
        }

        // Add to vector
        Position low(low_coords[0], low_coords[1], low_coords[2]);
        Position hi(hi_coords[0], hi_coords[1], hi_coords[2]);
        sources.push_back(std::make_shared<BoxSource>(low, hi, chi, weight));
      } else if (type == "point") {
        // Get coordinates
        std::vector<double> coords;
        if (input["sources"][s]["coords"] &&
            input["sources"][s]["coords"].IsSequence() &&
            input["sources"][s]["coords"].size() == 3) {
          coords = input["sources"][s]["coords"].as<std::vector<double>>();
        } else {
          std::string mssg = "No coordinates given for point source.";
          fatal_error(mssg, __FILE__, __LINE__);
        }

        // Add to vector
        Position point(coords[0], coords[1], coords[2]);
        sources.push_back(std::make_shared<PointSource>(point, chi, weight));
      } else {
        std::string mssg = "Unknown source type " + type + ".";
        fatal_error(mssg, __FILE__, __LINE__);
      }
    }
  } else if (!input["settings"]["insource"]) {
    // No source file given
    std::string mssg = "No source specified for problem.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
}

void make_simulation(YAML::Node input) {
  // Get simulation type
  std::string sim_type;
  if (input["settings"]["simulation"] &&
      input["settings"]["simulation"].IsScalar()) {
    sim_type = input["settings"]["simulation"].as<std::string>();
  } else {
    std::string mssg = "No simulation type provided.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  if (sim_type == "k-eigenvalue") {
    simulation = std::make_shared<PowerIterator>(tallies, settings, transporter,
                                                 sources);
  } else {
    std::string mssg = "Unknown simulation type " + sim_type + ".";
    fatal_error(mssg, __FILE__, __LINE__);
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