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
#include <cmath>
#include <geometry/hex_lattice.hpp>
#include <utils/error.hpp>

HexLattice::HexLattice(uint32_t nrings, uint32_t nz, double p, double pz,
                       double x, double y, double z, Top t, uint32_t i_id,
                       std::string i_name)
    : Lattice{i_id, i_name},
      pitch_{p},
      pitch_z_{pz},
      X_o{x},
      Y_o{y},
      Z_o{z},
      Nrings{nrings},
      Nz{nz},
      Nhex{0},
      width{0},
      mid_qr{0},
      top_{t} {
  Nhex = 0;
  for (size_t r = 0; r < Nrings; r++) {
    if (r == 0)
      Nhex += 1;
    else
      Nhex += 6 * r;
  }

  width = 2 * (Nrings - 1) + 1;

  mid_qr = width / 2;
}

bool HexLattice::is_inside(Position r, Direction /*u*/) const {
  // Get coordinates in frame of center tile
  Position r_o{r.x() - X_o, r.y() - Y_o, r.z() - Z_o};

  // Get qr of hearest hex tile
  std::array<int32_t, 2> qr = get_nearest_hex(r_o);

  // Get the ring of hex tile
  uint32_t ring = get_ring(qr);

  // See if valid ring or not
  if (ring >= Nrings)
    // Invalid ring
    return false;

  // Check z bin now
  double Z_low = Z_o - 0.5 * static_cast<double>(Nz) * pitch_z_;
  int32_t nz = std::floor((r.z() - Z_low) / pitch_z_);
  if (nz < 0 || nz >= static_cast<int32_t>(Nz)) return false;

  return true;
}

std::shared_ptr<Cell> HexLattice::get_cell(Position r, Direction u,
                                           int32_t on_surf) const {
  // Get coordinates in frame of center tile
  Position r_o{r.x() - X_o, r.y() - Y_o, r.z() - Z_o};

  // Get qr of hearest hex tile
  std::array<int32_t, 3> qrz = get_tile(r_o, u);

  // Get the ring of hex tile
  uint32_t ring = get_ring({qrz[0], qrz[1]});

  // See if valid ring or not
  if (ring >= Nrings) {
    // Invalid ring
    if (outer_universe_index == -1) {
      return nullptr;
    } else {
      return geometry::universes[outer_universe_index]->get_cell(r, u, on_surf);
    }
  }

  // Check z bin now
  if (qrz[2] < 0 || qrz[2] >= static_cast<int32_t>(Nz)) {
    if (outer_universe_index == -1) {
      return nullptr;
    } else {
      return geometry::universes[outer_universe_index]->get_cell(r, u, on_surf);
    }
  }

  // Inside a lattice bin
  size_t indx = linear_index({qrz[0], qrz[1]}, qrz[2]);

  // If -1, send to outside universe
  if (lattice_universes[indx] == -1) {
    if (outer_universe_index == -1) {
      return nullptr;
    } else {
      return geometry::universes[outer_universe_index]->get_cell(r, u, on_surf);
    }
  }

  // Move coordinates to tile center
  Position center = tile_center(qrz[0], qrz[1], qrz[2]);
  Position r_tile = r_o - center;
  return geometry::universes[lattice_universes[indx]]->get_cell(r_tile, u,
                                                                on_surf);
}

Cell* HexLattice::get_cell_naked_ptr(Position r, Direction u,
                                     int32_t on_surf) const {
  // Get coordinates in frame of center tile
  Position r_o{r.x() - X_o, r.y() - Y_o, r.z() - Z_o};

  // Get qr of hearest hex tile
  std::array<int32_t, 3> qrz = get_tile(r_o, u);

  // Get the ring of hex tile
  uint32_t ring = get_ring({qrz[0], qrz[1]});

  // See if valid ring or not
  if (ring >= Nrings) {
    // Invalid ring
    if (outer_universe_index == -1) {
      return nullptr;
    } else {
      return geometry::universes[outer_universe_index]->get_cell_naked_ptr(
          r, u, on_surf);
    }
  }

  // Check z bin now
  if (qrz[2] < 0 || qrz[2] >= static_cast<int32_t>(Nz)) {
    if (outer_universe_index == -1) {
      return nullptr;
    } else {
      return geometry::universes[outer_universe_index]->get_cell_naked_ptr(
          r, u, on_surf);
    }
  }

  // Inside a lattice bin
  size_t indx = linear_index({qrz[0], qrz[1]}, qrz[2]);

  // If -1, send to outside universe
  if (lattice_universes[indx] == -1) {
    if (outer_universe_index == -1) {
      return nullptr;
    } else {
      return geometry::universes[outer_universe_index]->get_cell_naked_ptr(
          r, u, on_surf);
    }
  }

  // Move coordinates to tile center
  Position center = tile_center(qrz[0], qrz[1], qrz[2]);
  Position r_tile = r_o - center;
  return geometry::universes[lattice_universes[indx]]->get_cell_naked_ptr(
      r_tile, u, on_surf);
}

std::shared_ptr<Cell> HexLattice::get_cell(std::vector<GeoLilyPad>& stack,
                                           Position r, Direction u,
                                           int32_t on_surf) const {
  // Get coordinates in frame of center tile
  Position r_o{r.x() - X_o, r.y() - Y_o, r.z() - Z_o};

  // Get qr of hearest hex tile
  std::array<int32_t, 3> qrz = get_tile(r_o, u);

  // Get the ring of hex tile
  uint32_t ring = get_ring({qrz[0], qrz[1]});

  // See if valid ring or not
  if (ring >= Nrings) {
    // Invalid ring
    if (outer_universe_index == -1) {
      // Save info to stack
      stack.push_back({GeoLilyPad::PadType::Lattice, id_, r, qrz, false});
      return nullptr;
    } else {
      // Save info to stack
      stack.push_back({GeoLilyPad::PadType::Lattice, id_, r, qrz, true});
      return geometry::universes[outer_universe_index]->get_cell(stack, r, u,
                                                                 on_surf);
    }
  }

  // Check z bin now
  if (qrz[2] < 0 || qrz[2] >= static_cast<int32_t>(Nz)) {
    if (outer_universe_index == -1) {
      // Save info to stack
      stack.push_back({GeoLilyPad::PadType::Lattice, id_, r, qrz, false});
      return nullptr;
    } else {
      // Save info to stack
      stack.push_back({GeoLilyPad::PadType::Lattice, id_, r, qrz, true});
      return geometry::universes[outer_universe_index]->get_cell(stack, r, u,
                                                                 on_surf);
    }
  }

  // Inside a lattice bin
  size_t indx = linear_index({qrz[0], qrz[1]}, qrz[2]);

  // If -1, send to outside universe
  if (lattice_universes[indx] == -1) {
    if (outer_universe_index == -1) {
      // Save info to stack
      stack.push_back({GeoLilyPad::PadType::Lattice, id_, r, qrz, false});
      return nullptr;
    } else {
      // Save info to stack
      stack.push_back({GeoLilyPad::PadType::Lattice, id_, r, qrz, true});
      return geometry::universes[outer_universe_index]->get_cell(stack, r, u,
                                                                 on_surf);
    }
  }

  // Move coordinates to tile center
  Position center = tile_center(qrz[0], qrz[1], qrz[2]);
  Position r_tile = r_o - center;
  // Save info to stack
  stack.push_back({GeoLilyPad::PadType::Lattice, id_, r, qrz, false});
  return geometry::universes[lattice_universes[indx]]->get_cell(stack, r_tile,
                                                                u, on_surf);
}

void HexLattice::set_elements(std::vector<int32_t> univs) {
  // Make sure proper number of universes are provided
  if (univs.size() != Nhex * Nz) {
    std::string mssg = "Improper number of universes for HexLattice.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  lattice_universes.resize(width * width * Nz, -1);

  size_t indx = 0;
  for (size_t az = 0; az < Nz; az++) {
    for (size_t ar = 0; ar < width; ar++) {
      for (size_t aq = 0; aq < width; aq++) {
        int32_t q = aq - mid_qr;
        int32_t r = ar - mid_qr;

        uint32_t ring = get_ring({q, r});

        if (ring < Nrings) {
          lattice_universes[linear_index({q, r}, az)] = univs[indx];
          indx += 1;
        }
      }
    }
  }
}

size_t HexLattice::linear_index(std::array<int32_t, 2> qr, uint32_t z) const {
  // In qr coordinates shifted by mid_qr, both q and r can range from 0 to
  // width - 1. This then forms a width*width square array.
  size_t indx = 0;
  uint32_t q = qr[0] + mid_qr;
  uint32_t r = qr[1] + mid_qr;
  indx = z * (width * width) + r * (width) + q;

  if (indx >= lattice_universes.size()) {
    std::string mssg = "Invalid index for HexLattice.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  return indx;
}

std::array<int32_t, 2> HexLattice::get_nearest_hex(Position p) const {
  double q, r, det;
  if (top_ == Top::Flat) {
    det = -pitch_ * pitch_ * sin_pi_3;

    q = (pitch_ / det) * (0. * p.x() - 1. * p.y());
    r = (pitch_ / det) * (-sin_pi_3 * p.x() + cos_pi_3 * p.y());
  } else {
    det = -pitch_ * pitch_ * cos_pi_6;

    q = (pitch_ / det) * (sin_pi_6 * p.x() - cos_pi_6 * p.y());
    r = (pitch_ / det) * (-1. * p.x() - 0. * p.y());
  }

  // Do rounding magic to get nearest hexagon
  // First convert to cube coordinates
  double x = q;
  double z = r;
  double y = -x - z;

  // Do rounding
  double rx = std::round(x);
  double ry = std::round(y);
  double rz = std::round(z);

  // Differences
  double x_diff = std::abs(rx - x);
  double y_diff = std::abs(ry - y);
  double z_diff = std::abs(rz - z);

  // Check three cases, correct rounded coordinates
  if (x_diff > y_diff && x_diff > z_diff) {
    rx = -ry - rz;
  } else if (y_diff > x_diff && y_diff > z_diff) {
    ry = -rx - rz;
  } else {
    rz = -rx - ry;
  }

  // Convert back to axial coordinates for return
  return {static_cast<int32_t>(rx), static_cast<int32_t>(rz)};
}

std::array<int32_t, 3> HexLattice::get_tile(Position p, Direction /*u*/) const {
  std::array<int32_t, 2> qr = get_nearest_hex(p);
  double Z_low = Z_o - 0.5 * static_cast<double>(Nz) * pitch_z_;
  int32_t nz = std::floor((p.z() - Z_low) / pitch_z_);
  return {qr[0], qr[1], nz};
}

Position HexLattice::get_hex_center(std::array<int32_t, 2> qr) const {
  double x, y;
  double q = static_cast<double>(qr[0]);
  double r = static_cast<double>(qr[1]);

  if (top_ == Top::Flat) {
    x = pitch_ * (cos_pi_3 * q + 1. * r);
    y = pitch_ * (sin_pi_3 * q + 0. * r);
  } else {
    x = pitch_ * (0. * q + cos_pi_6 * r);
    y = pitch_ * (1. * q + sin_pi_6 * r);
  }

  return Position(x, y, 0.);
}

Position HexLattice::tile_center(int q, int r, int nz) const {
  double Z_low = Z_o - 0.5 * static_cast<double>(Nz) * pitch_z_;
  Position r_hex = get_hex_center({q, r});
  double z = (static_cast<double>(nz) + 0.5) * pitch_z_ + Z_low;
  return Position(r_hex.x(), r_hex.y(), z);
}

uint32_t HexLattice::get_ring(std::array<int32_t, 2> qr) const {
  int32_t x = qr[0];
  int32_t z = qr[1];
  int32_t y = -x - z;

  uint32_t ax = static_cast<uint32_t>(std::abs(x));
  uint32_t ay = static_cast<uint32_t>(std::abs(y));
  uint32_t az = static_cast<uint32_t>(std::abs(z));

  // Ring is the max abs of x, y z
  uint32_t temp_abs = std::max(ax, ay);
  return std::max(temp_abs, az);
}

double HexLattice::distance_to_tile_boundary(
    Position r_local, Direction u, std::array<int32_t, 3> tile) const {
  Position center = tile_center(tile[0], tile[1], tile[2]);
  Position r_tile = r_local - center;

  if (top_ == Top::Pointy)
    return distance_to_tile_boundary_pointy(r_tile, u);
  else
    return distance_to_tile_boundary_flat(r_tile, u);
}

double HexLattice::distance_to_tile_boundary_pointy(Position r_tile,
                                                    Direction u) const {
  // Points needed to calculate slopes for lines. Only need two sets
  // because of symmetry reasons.
  //    (x1,y1)
  /*      /\
         |  |(x2,y2)  */
  double x1 = 0.;
  double y1 = pitch_ / (2. * cos_pi_6);
  double x2 = pitch_ / 2.;
  double y2 = pitch_ / (2. * sin_pi_6);

  double d = INF;

  // Get distance to all six hexagon surfaces
  double d1 = distance_to_line(r_tile, u, x1, y1, x2, y1);
  double d2 = distance_to_line(r_tile, u, x2, y2, x2, -y2);
  double d3 = distance_to_line(r_tile, u, x2, -y2, x1, -y1);
  double d4 = distance_to_line(r_tile, u, x1, -y1, -x2, -y2);
  double d5 = distance_to_line(r_tile, u, -x2, -y2, -x2, y2);
  double d6 = distance_to_line(r_tile, u, -x2, y2, x1, y1);

  // Get distance to z planes
  double dzl = (-pitch_z_ * 0.5 - r_tile.z()) / u.z();
  double dzu = (pitch_z_ * 0.5 - r_tile.z()) / u.z();

  if (d1 > 0. && d1 < d) d = d1;
  if (d2 > 0. && d2 < d) d = d2;
  if (d3 > 0. && d3 < d) d = d3;
  if (d4 > 0. && d4 < d) d = d4;
  if (d5 > 0. && d5 < d) d = d5;
  if (d6 > 0. && d6 < d) d = d6;
  if (dzl > 0. && dzl < d) d = dzl;
  if (dzu > 0. && dzu < d) d = dzu;

  return d;
}

double HexLattice::distance_to_tile_boundary_flat(Position r_tile,
                                                  Direction u) const {
  // Points needed to calculate slopes for lines. Only need two sets
  // because of symmetry reasons.
  //
  /*        --- (x1, y1)
           /   \
                 (x2, y2)  */
  double x1 = pitch_ / (2. * sin_pi_6);
  double y1 = pitch_ / 2.;
  double x2 = pitch_ / (2. * cos_pi_6);
  double y2 = 0.;

  double d = INF;

  // Get distance to all six hexagon surfaces
  double d1 = distance_to_line(r_tile, u, x1, y1, x2, y1);
  double d2 = distance_to_line(r_tile, u, x2, y2, x1, -y1);
  double d3 = distance_to_line(r_tile, u, x1, -y1, -x1, -y1);
  double d4 = distance_to_line(r_tile, u, -x1, -y1, -x2, -y2);
  double d5 = distance_to_line(r_tile, u, -x2, -y2, -x1, y1);
  double d6 = distance_to_line(r_tile, u, -x1, y1, x1, y1);

  // Get distance to z planes
  double dzl = (-pitch_z_ * 0.5 - r_tile.z()) / u.z();
  double dzu = (pitch_z_ * 0.5 - r_tile.z()) / u.z();

  if (d1 > 0. && d1 < d) d = d1;
  if (d2 > 0. && d2 < d) d = d2;
  if (d3 > 0. && d3 < d) d = d3;
  if (d4 > 0. && d4 < d) d = d4;
  if (d5 > 0. && d5 < d) d = d5;
  if (d6 > 0. && d6 < d) d = d6;
  if (dzl > 0. && dzl < d) d = dzl;
  if (dzu > 0. && dzu < d) d = dzu;

  return d;
}

double HexLattice::distance_to_line(Position r, Direction u, double x1,
                                    double y1, double x2, double y2) const {
  // Get A, B, and D for plane from y - y1 = m (x - x1) ; m = (y2 - y1)/(x2 -
  // x1)
  double A = y2 - y1;
  double B = x1 - x2;
  double D = (x2 - x1) * y1 - (y2 - y1) * x1;

  double num = D - A * r.x() - B * r.y();
  double denom = A * u.x() + B * u.y();

  double d = num / denom;

  if (d < 0.)
    return INFINITY;
  else
    return d;
}

void make_hex_lattice(YAML::Node latt_node, YAML::Node input) {
  // Get id
  uint32_t id;
  if (latt_node["id"] && latt_node["id"].IsScalar()) {
    id = latt_node["id"].as<uint32_t>();
  } else {
    std::string mssg = "Lattice must have a valid id.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get name if present
  std::string name = "";
  if (latt_node["name"] && latt_node["name"].IsScalar()) {
    name = latt_node["name"].as<std::string>();
  }

  // Get shape
  std::vector<uint32_t> shape;
  shape.resize(2);
  if (latt_node["shape"] && latt_node["shape"].IsSequence() &&
      latt_node["shape"].size() == 2) {
    for (size_t s = 0; s < 2; s++) {
      shape[s] = latt_node["shape"][s].as<uint32_t>();
    }
  } else {
    std::string mssg = "Lattice must have a valid shape.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get pitch
  std::vector<double> pitch;
  pitch.resize(2);
  if (latt_node["pitch"] && latt_node["pitch"].IsSequence() &&
      latt_node["pitch"].size() == 2) {
    for (size_t s = 0; s < 2; s++) {
      pitch[s] = latt_node["pitch"][s].as<double>();
    }
  } else {
    std::string mssg = "Lattice must have a valid pitch.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get origin
  std::vector<double> origin;
  origin.resize(3);
  if (latt_node["origin"] && latt_node["origin"].IsSequence() &&
      latt_node["origin"].size() == 3) {
    for (size_t s = 0; s < 3; s++) {
      origin[s] = latt_node["origin"][s].as<double>();
    }
  } else {
    std::string mssg = "Lattice must have a valid origin.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Vector for universes
  std::vector<int32_t> uni_indicies;
  if (latt_node["universes"] && latt_node["universes"].IsSequence()) {
    // Go through and check each universe
    for (size_t u = 0; u < latt_node["universes"].size(); u++) {
      int32_t u_id = latt_node["universes"][u].as<int32_t>();
      if (u_id == -1) {
        uni_indicies.push_back(u_id);
      } else if (universe_id_to_indx.find(u_id) == universe_id_to_indx.end()) {
        // Need to find universe
        find_universe(input, u_id);
        uni_indicies.push_back(universe_id_to_indx[u_id]);
      } else {
        uni_indicies.push_back(universe_id_to_indx[u_id]);
      }
    }
  } else {
    std::string mssg =
        "Lattice instance must have a valid universes definition.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make sure lattice id not taken
  if (lattice_id_to_indx.find(id) != lattice_id_to_indx.end()) {
    std::string mssg = "Lattice id " + std::to_string(id) +
                       " appears"
                       " multiple times.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get top type
  HexLattice::Top top = HexLattice::Top::Pointy;
  if (latt_node["top"] && latt_node["top"].IsScalar()) {
    if (latt_node["top"].as<std::string>() == "pointy") {
      top = HexLattice::Top::Pointy;
    } else if (latt_node["top"].as<std::string>() == "flat") {
      top = HexLattice::Top::Flat;
    } else {
      std::string mssg = " Uknown top for hexagonal lattice " +
                         latt_node["top"].as<std::string>() + ".";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  }

  // Make lattice
  std::shared_ptr<Lattice> lat = std::make_shared<HexLattice>(
      shape[0], shape[1], pitch[0], pitch[1], origin[0], origin[1], origin[2],
      top, id, name);

  // Set universes
  lat->set_elements(uni_indicies);

  // Get outside
  int32_t out_id;
  if (latt_node["outer"] && latt_node["outer"].IsScalar()) {
    out_id = latt_node["outer"].as<int32_t>();
    if (out_id != -1) {
      // Find outside universe
      if (universe_id_to_indx.find(out_id) == universe_id_to_indx.end()) {
        // Need to find universe
        find_universe(input, out_id);
      }
      lat->set_outisde_universe(universe_id_to_indx[out_id]);
    }
  }

  // Set lattice
  lattice_id_to_indx[id] = geometry::lattices.size();
  geometry::lattices.push_back(lat);
}
