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
#include <simulation/collision_mesh_tally.hpp>
#include <simulation/mesh_tally.hpp>
#include <simulation/track_length_mesh_tally.hpp>
#include <utils/error.hpp>
#include <utils/mpi.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>
#include <vector>

MeshTally::MeshTally(Position low, Position hi, uint64_t nx, uint64_t ny,
                     uint64_t nz, const std::vector<double> &ebounds,
                     TallyQuantity q, std::string fname, uint32_t mt)
    : r_low{low},
      r_hi{hi},
      Nx{nx},
      Ny{ny},
      Nz{nz},
      dx(),
      dy(),
      dz(),
      g(),
      net_weight(1.),
      energy_bounds(ebounds),
      quantity(q),
      fname(fname),
      mt(mt),
      tally_gen(),
      tally_avg(),
      tally_var() {
  dx = (r_hi.x() - r_low.x()) / static_cast<double>(Nx);
  dy = (r_hi.y() - r_low.y()) / static_cast<double>(Ny);
  dz = (r_hi.z() - r_low.z()) / static_cast<double>(Nz);

  uint32_t Ne = energy_bounds.size() - 1;

  // Allocate and fill arrays to zero
  tally_gen.reallocate({Ne, Nx, Ny, Nz});
  tally_gen.fill(0.);

  // Only allocate average and variance if we are the master !
  if (mpi::rank == 0) {
    tally_avg.reallocate({Ne, Nx, Ny, Nz});
    tally_avg.fill(0.);

    tally_var.reallocate({Ne, Nx, Ny, Nz});
    tally_var.fill(0.);
  }
}

void MeshTally::set_net_weight(double W) { net_weight = W; }

void MeshTally::record_generation(double gen, double multiplier) {
  g = gen;

  // All worker threads must send their generation score to the master.
  // Master must recieve all generations scores from workers and add
  // them to it's own generation score.
  mpi::Reduce_sum(&tally_gen[0], static_cast<int>(tally_gen.size()), 0);

  // Only try to update average and variance is we are master, as worker
  // processes don't have copies of this data, so it will seg-fault.
  if (mpi::rank == 0) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (size_t i = 0; i < tally_gen.size(); i++) {
      // Get new average
      double old_avg = tally_avg[i];
      double val = tally_gen[i] * multiplier;
      double avg = old_avg + (val - old_avg) / gen;
      tally_avg[i] = avg;

      // Get new variance
      double var = tally_var[i];
      var = var + (((val - old_avg) * (val - avg) - (var)) / gen);
      tally_var[i] = var;
    }
  }
}

void MeshTally::clear_generation() { tally_gen.fill(0.); }

void MeshTally::write_tally() {
  // Only master can write tallies, as only master has a copy
  // of the mean and variance.
  if (mpi::rank == 0) {
    std::string mssg = " Writing " + fname + " tally file...\n";
    Output::instance()->write(mssg);

    std::ofstream file(fname + "_mesh.txt");

    // First write coordinates and number of groups
    file << " X:";
    for (int i = 0; i < static_cast<int>(Nx) - 1; i++) {
      double x = i * dx + (0.5 * dx) + r_low.x();
      file << x << ",";
    }
    file << (Nx - 1) * dx + (0.5 * dx) + r_low.x() << "\n";

    file << " Y:";
    for (int i = 0; i < static_cast<int>(Ny) - 1; i++) {
      double y = i * dy + (0.5 * dy) + r_low.y();
      file << y << ",";
    }
    file << (Ny - 1) * dy + (0.5 * dy) + r_low.y() << "\n";

    file << " Z:";
    for (int i = 0; i < static_cast<int>(Nz) - 1; i++) {
      double z = i * dz + (0.5 * dz) + r_low.z();
      file << z << ",";
    }
    file << (Nz - 1) * dz + (0.5 * dz) + r_low.z() << "\n";

    file << " ENERGY:";
    for (size_t i = 0; i < energy_bounds.size() - 1; i++) {
      file << energy_bounds[i] << ",";
    }
    file << energy_bounds.back() << "\n";

    file << " QUANTITY: ";

    switch (quantity) {
      case TallyQuantity::Flux:
        file << "Flux\n";
        break;

      case TallyQuantity::Elastic:
        file << "Elastic\n";
        break;

      case TallyQuantity::Absorption:
        file << "Absorption\n";
        break;

      case TallyQuantity::Fission:
        file << "Fission\n";
        break;

      case TallyQuantity::Total:
        file << "Total\n";
        break;

      case TallyQuantity::MT:
        file << "MT = " << mt << "\n";
        break;

      case TallyQuantity::RealFlux:
        file << "RealFlux\n";
        break;

      case TallyQuantity::RealElastic:
        file << "RealElastic\n";
        break;

      case TallyQuantity::RealAbsorption:
        file << "RealAbsorption\n";
        break;

      case TallyQuantity::RealFission:
        file << "RealFission\n";
        break;

      case TallyQuantity::RealTotal:
        file << "RealTotal\n";
        break;

      case TallyQuantity::RealMT:
        file << "RealMT = " << mt << "\n";
        break;

      case TallyQuantity::ImgFlux:
        file << "ImgFlux\n";
        break;

      case TallyQuantity::ImgElastic:
        file << "ImgElastic\n";
        break;

      case TallyQuantity::ImgAbsorption:
        file << "ImgAbsorption\n";
        break;

      case TallyQuantity::ImgFission:
        file << "ImgFission\n";
        break;

      case TallyQuantity::ImgTotal:
        file << "ImgTotal\n";
        break;

      case TallyQuantity::ImgMT:
        file << "ImgMT = " << mt << "\n";
        break;

      case TallyQuantity::MagFlux:
        file << "MagFlux\n";
        break;

      case TallyQuantity::MagSqrFlux:
        file << "MagSqrFlux\n";
        break;
    }

    file.close();

    // Convert flux_var to the error on the mean
    for (size_t l = 0; l < tally_var.size(); l++)
      tally_var[l] = std::sqrt(tally_var[l] / g);

    tally_avg.save(fname + "_avg.npy");
    tally_var.save(fname + "_err.npy");
  }
}

std::shared_ptr<MeshTally> make_mesh_tally(const YAML::Node &node) {
  // First get type of estimator. Default is collision
  std::string estimator_str = "collision";
  if (node["estimator"]) {
    estimator_str = node["estimator"].as<std::string>();
  }

  // Get low position
  double xl = 0., yl = 0., zl = 0.;
  if (!node["low"] || !node["low"].IsSequence() || node["low"].size() != 3) {
    std::string mssg = "Now valid low entry for mesh tally.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  xl = node["low"][0].as<double>();
  yl = node["low"][1].as<double>();
  zl = node["low"][2].as<double>();
  Position plow(xl, yl, zl);

  // Get hi position
  double xh = 0., yh = 0., zh = 0.;
  if (!node["hi"] || !node["hi"].IsSequence() || node["hi"].size() != 3) {
    std::string mssg = "Now valid hi entry for mesh tally.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  xh = node["hi"][0].as<double>();
  yh = node["hi"][1].as<double>();
  zh = node["hi"][2].as<double>();
  Position phi(xh, yh, zh);

  // Check positions
  if (plow.x() >= phi.x() || plow.y() >= phi.y() || plow.x() >= phi.x()) {
    std::string mssg =
        "Low coordinates for mesh tally must be less than hi coordinates.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get shape
  uint64_t nx = 1, ny = 1, nz = 1;
  if (node["shape"] &&
      (!node["shape"].IsSequence() || node["shape"].size() != 3)) {
    std::string mssg = "Invalid shape provided to mesh tally.";
  } else if (node["shape"]) {
    int64_t tmp_nx = node["shape"][0].as<int64_t>();
    int64_t tmp_ny = node["shape"][1].as<int64_t>();
    int64_t tmp_nz = node["shape"][2].as<int64_t>();

    if (tmp_nx < 1 || tmp_ny < 1 || tmp_nz < 1) {
      std::string mssg =
          "Mesh tally shapes must be values greater than or equal to 1.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    nx = static_cast<uint64_t>(tmp_nx);
    ny = static_cast<uint64_t>(tmp_ny);
    nz = static_cast<uint64_t>(tmp_nz);
  }

  // Get energy bounds
  std::vector<double> ebounds;
  if (!node["energy-bounds"] || !node["energy-bounds"].IsSequence()) {
    std::string mssg = "No valid energy-bounds enetry provided to mesh tally.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  ebounds = node["energy-bounds"].as<std::vector<double>>();
  if (ebounds.size() < 2) {
    std::string mssg = "Energy-bounds must have at least two entries.";
    fatal_error(mssg, __FILE__, __LINE__);
  } else if (!std::is_sorted(ebounds.begin(), ebounds.end())) {
    std::string mssg = "Energy-bounds must be sorted.";
    fatal_error(mssg, __FILE__, __LINE__);
  } else if (ebounds.front() < 0.) {
    std::string mssg = "All energy-bounds entries must be positive.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get name
  std::string fname;
  if (!node["name"] || !node["name"].IsScalar()) {
    std::string mssg = "No valid name provided to mesh tally.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  fname = node["name"].as<std::string>();

  // Get tally quantity
  uint32_t mt = 0;
  std::string quant_str = "none";
  MeshTally::TallyQuantity quantity = MeshTally::TallyQuantity::Flux;
  if (!node["quantity"] || !node["quantity"].IsScalar()) {
    std::string mssg = "No quantity entry provided to mesh tally.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  quant_str = node["quantity"].as<std::string>();
  if (quant_str == "flux") {
    quantity = MeshTally::TallyQuantity::Flux;
  } else if (quant_str == "total") {
    quantity = MeshTally::TallyQuantity::Total;
  } else if (quant_str == "elastic") {
    quantity = MeshTally::TallyQuantity::Elastic;
  } else if (quant_str == "absorption") {
    quantity = MeshTally::TallyQuantity::Absorption;
  } else if (quant_str == "fission") {
    quantity = MeshTally::TallyQuantity::Fission;
  } else if (quant_str == "mt") {
    quantity = MeshTally::TallyQuantity::MT;

    if (settings::energy_mode == settings::EnergyMode::MG) {
      // Can't do an MT tally in MG mode !
      std::string mssg = "Cannot do MT tallies in multi-group mode.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Check for mt
    if (!node["mt"] || !node["mt"].IsScalar()) {
      std::string mssg =
          "Quantity of \"mt\" selected, but no provided mt value.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    int32_t tmp_mt = node["mt"].as<int32_t>();
    if (tmp_mt < 4 || tmp_mt > 891) {
      std::string mssg =
          "The value " + std::to_string(tmp_mt) + " is not a valid MT.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    mt = static_cast<uint32_t>(tmp_mt);
  } else if (quant_str == "real-flux") {
    quantity = MeshTally::TallyQuantity::RealFlux;
  } else if (quant_str == "real-total") {
    quantity = MeshTally::TallyQuantity::RealTotal;
  } else if (quant_str == "real-elastic") {
    quantity = MeshTally::TallyQuantity::RealElastic;
  } else if (quant_str == "real-absorption") {
    quantity = MeshTally::TallyQuantity::RealAbsorption;
  } else if (quant_str == "real-fission") {
    quantity = MeshTally::TallyQuantity::RealFission;
  } else if (quant_str == "real-mt") {
    quantity = MeshTally::TallyQuantity::RealMT;

    if (settings::energy_mode == settings::EnergyMode::MG) {
      // Can't do an MT tally in MG mode !
      std::string mssg = "Cannot do MT tallies in multi-group mode.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Check for mt
    if (!node["mt"] || !node["mt"].IsScalar()) {
      std::string mssg =
          "Quantity of \"mt\" selected, but no provided mt value.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    int32_t tmp_mt = node["mt"].as<int32_t>();
    if (tmp_mt < 4 || tmp_mt > 891) {
      std::string mssg =
          "The value " + std::to_string(tmp_mt) + " is not a valid MT.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    mt = static_cast<uint32_t>(tmp_mt);
  } else if (quant_str == "img-flux") {
    quantity = MeshTally::TallyQuantity::ImgFlux;
  } else if (quant_str == "img-total") {
    quantity = MeshTally::TallyQuantity::ImgTotal;
  } else if (quant_str == "img-elastic") {
    quantity = MeshTally::TallyQuantity::ImgElastic;
  } else if (quant_str == "img-absorption") {
    quantity = MeshTally::TallyQuantity::ImgAbsorption;
  } else if (quant_str == "img-fission") {
    quantity = MeshTally::TallyQuantity::ImgFission;
  } else if (quant_str == "img-mt") {
    quantity = MeshTally::TallyQuantity::ImgMT;

    if (settings::energy_mode == settings::EnergyMode::MG) {
      // Can't do an MT tally in MG mode !
      std::string mssg = "Cannot do MT tallies in multi-group mode.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Check for mt
    if (!node["mt"] || !node["mt"].IsScalar()) {
      std::string mssg =
          "Quantity of \"mt\" selected, but no provided mt value.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    int32_t tmp_mt = node["mt"].as<int32_t>();
    if (tmp_mt < 4 || tmp_mt > 891) {
      std::string mssg =
          "The value " + std::to_string(tmp_mt) + " is not a valid MT.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    mt = static_cast<uint32_t>(tmp_mt);
  } else if (quant_str == "mag-flux") {
    quantity = MeshTally::TallyQuantity::MagFlux;
  } else if (quant_str == "mag-sqr-flux") {
    quantity = MeshTally::TallyQuantity::MagSqrFlux;
  } else {
    std::string mssg = "Unkown tally quantity \"" + quant_str + "\".";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Construct based on estimator type
  if (estimator_str == "collision") {
    return std::make_shared<CollisionMeshTally>(plow, phi, nx, ny, nz, ebounds,
                                                quantity, fname, mt);
  } else if (estimator_str == "track-length") {
    // TODO send settings along too, and throw error if trying to use TLE
    // with DT or CT
    if ((settings::tracking == settings::TrackingMode::DELTA_TRACKING ||
         settings::tracking == settings::TrackingMode::CARTER_TRACKING) &&
        (quantity != MeshTally::TallyQuantity::Flux &&
         quantity != MeshTally::TallyQuantity::RealFlux &&
         quantity != MeshTally::TallyQuantity::ImgFlux &&
         quantity != MeshTally::TallyQuantity::MagFlux &&
         quantity != MeshTally::TallyQuantity::MagSqrFlux)) {
      std::string mssg =
          "Cannot use track-length estimators for non-flux quantities with "
          "delta-tracking or carter-tracking.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    return std::make_shared<TrackLengthMeshTally>(plow, phi, nx, ny, nz,
                                                  ebounds, quantity, fname, mt);
    return nullptr;
  } else {
    std::string mssg = "Unknown estimator type of \"" + estimator_str + "\".";
    fatal_error(mssg, __FILE__, __LINE__);
    return nullptr;
  }
}
