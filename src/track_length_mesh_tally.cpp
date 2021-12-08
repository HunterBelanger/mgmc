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
#include <simulation/track_length_mesh_tally.hpp>
#include <utils/constants.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>

void TrackLengthMeshTally::score_flight(const Particle &p, double d,
                                        MaterialHelper &mat) {
  // Local position and direction copies which we will use for tallying
  Position r = p.r();
  Direction u = p.u();

  // Initialize the spatial indices with the starting position
  // of the particle
  int i = 0, j = 0, k = 0;
  std::array<int, 3> on{0, 0, 0};
  initialize_indices(r, u, i, j, k, on);

  // First, check if we are inside the tally region
  bool inside_tally_region = false;
  if (i >= 0 && i < static_cast<int>(Nx) && j >= 0 &&
      j < static_cast<int>(Ny) && k >= 0 && k < static_cast<int>(Nz)) {
    inside_tally_region = true;
  }

  // If we aren't inside, are we going to hit it along the trajectory ?
  // If so, we can fast-forward to that position.
  if (!inside_tally_region) {
    if (!find_entry_point(r, u, d)) {
      return;
    }
    inside_tally_region = true;
  }

  /*if (p.wgt() == 0. && p.wgt2() == 0.) {
    std::stringstream out;
    out << " Zero weights in TL.\n";
    std::cout << out.str();
  }*/

  // Calculate base score, absed on the quantity
  double base_score = 1. / net_weight;
  // Must multiply score by correct factor depending on the observed
  // quantity. q = 0 corresponds to just the flux, so no modification.
  switch (quantity) {
    case TallyQuantity::Flux:
      base_score *= p.wgt();
      break;

    case TallyQuantity::Elastic:
      base_score *= p.wgt() * mat.Eelastic(p.E());
      break;

    case TallyQuantity::Absorption:
      base_score *= p.wgt() * mat.Ea(p.E());
      break;

    case TallyQuantity::Fission:
      base_score *= p.wgt() * mat.Ef(p.E());
      break;

    case TallyQuantity::Total:
      base_score *= p.wgt() * mat.Et(p.E());
      break;

    case TallyQuantity::MT:
      base_score *= p.wgt() * mat.Emt(mt, p.E());
      break;

    case TallyQuantity::RealFlux:
      base_score *= p.wgt();
      break;

    case TallyQuantity::RealElastic:
      base_score *= p.wgt() * mat.Eelastic(p.E());
      break;

    case TallyQuantity::RealAbsorption:
      base_score *= p.wgt() * mat.Ea(p.E());
      break;

    case TallyQuantity::RealFission:
      base_score *= p.wgt() * mat.Ef(p.E());
      break;

    case TallyQuantity::RealTotal:
      base_score *= p.wgt() * mat.Et(p.E());
      break;

    case TallyQuantity::RealMT:
      base_score *= p.wgt() * mat.Emt(mt, p.E());
      break;

    case TallyQuantity::ImgFlux:
      base_score *= p.wgt2();
      break;

    case TallyQuantity::ImgElastic:
      base_score *= p.wgt2() * mat.Eelastic(p.E());
      break;

    case TallyQuantity::ImgAbsorption:
      base_score *= p.wgt2() * mat.Ea(p.E());
      break;

    case TallyQuantity::ImgFission:
      base_score *= p.wgt2() * mat.Ef(p.E());
      break;

    case TallyQuantity::ImgTotal:
      base_score *= p.wgt2() * mat.Et(p.E());
      break;

    case TallyQuantity::ImgMT:
      base_score *= p.wgt2() * mat.Emt(mt, p.E());
      break;

    case TallyQuantity::MagFlux:
      base_score *= std::sqrt(p.wgt() * p.wgt() + p.wgt2() * p.wgt2());
      break;

    case TallyQuantity::MagSqrFlux:
      base_score *= (p.wgt() * p.wgt() + p.wgt2() * p.wgt2());
      break;
  }

  // Get energy index with linear search
  int l = -1;
  for (size_t e = 0; e < energy_bounds.size() - 1; e++) {
    if (energy_bounds[e] <= p.E() && p.E() <= energy_bounds[e + 1]) {
      l = static_cast<int>(e);
      break;
    }
  }
  // Don't score anything if we didn't find an energy bin
  if (l == -1) {
    return;
  }
  uint64_t uE = static_cast<uint64_t>(l);

  // Distance remaining to tally
  double distance_remaining = d;

  while (true) {
    // Distance we will travel in this cell
    auto next_tile = distance_to_next_index(r, u, i, j, k, on);

    if (next_tile.first == INF) {
      // Something went wrong.... Don't score.
      Output::instance()->save_warning("Problem encountered with mesh tally " +
                                       fname + ".");
      break;
    }

    double d_tile = std::min(next_tile.first, distance_remaining);

    // Make the score if we are in a valid cell
    if (i >= 0 && i < static_cast<int>(Nx) && j >= 0 &&
        j < static_cast<int>(Ny) && k >= 0 && k < static_cast<int>(Nz)) {
      uint64_t ui = static_cast<uint64_t>(i);
      uint64_t uj = static_cast<uint64_t>(j);
      uint64_t uk = static_cast<uint64_t>(k);
#ifdef _OPENMP
#pragma omp atomic
#endif
      tally_gen(uE, ui, uj, uk) += d_tile * base_score;
    } else {
      // If we arrive here, it means that we have left the tally region
      // when were we initially inside it. We can return here, as it's
      // impossible to go back in.
      return;
    }

    // Remove the traveled distance
    distance_remaining -= d_tile;

    if (distance_remaining <= 0.) break;

    // Update the position and cell indices
    r = r + d_tile * u;
    update_indices(next_tile.second, i, j, k, on);

  }  // While we still have to travel
}

bool TrackLengthMeshTally::find_entry_point(Position &r, const Direction &u,
                                            double &d_flight) const {
  // Set the initial distance to entry, which is going to be INF
  double dist = INF;

  // Get the distance to all bounding surfaces
  double d_xl = (r_low.x() - r.x()) / u.x();
  double d_yl = (r_low.y() - r.y()) / u.y();
  double d_zl = (r_low.z() - r.z()) / u.z();
  double d_xh = (r_hi.x() - r.x()) / u.x();
  double d_yh = (r_hi.y() - r.y()) / u.y();
  double d_zh = (r_hi.z() - r.z()) / u.z();

  if (d_xl > 0. && d_xl < dist && d_xl > 100 * SURFACE_COINCIDENT) dist = d_xl;
  if (d_xh > 0. && d_xh < dist && d_xh > 100 * SURFACE_COINCIDENT) dist = d_xh;
  if (d_yl > 0. && d_yl < dist && d_yl > 100 * SURFACE_COINCIDENT) dist = d_yl;
  if (d_yh > 0. && d_yh < dist && d_yh > 100 * SURFACE_COINCIDENT) dist = d_yh;
  if (d_zl > 0. && d_zl < dist && d_zl > 100 * SURFACE_COINCIDENT) dist = d_zl;
  if (d_zh > 0. && d_zh < dist && d_zh > 100 * SURFACE_COINCIDENT) dist = d_zh;

  if (dist != INF && dist < d_flight) {
    r = r + dist * u;
    d_flight -= dist;
    return true;
  }

  return false;
}

void TrackLengthMeshTally::initialize_indices(const Position &r,
                                              const Direction &u, int &i,
                                              int &j, int &k,
                                              std::array<int, 3> &on) {
  i = std::floor((r.x() - r_low.x()) / dx);
  j = std::floor((r.y() - r_low.y()) / dy);
  k = std::floor((r.z() - r_low.z()) / dz);
  on.fill(0);

  // Must handle case of being on a tile boundary
  // Get position at center of current tile
  double xc = r_low.x() + i * dx + 0.5 * dx;
  double yc = r_low.y() + j * dy + 0.5 * dy;
  double zc = r_low.z() + k * dz + 0.5 * dz;

  // Get tile boundaries
  double xl = xc - 0.5 * dx;
  double xh = xc + 0.5 * dx;
  double yl = yc - 0.5 * dy;
  double yh = yc + 0.5 * dy;
  double zl = zc - 0.5 * dz;
  double zh = zc + 0.5 * dz;

  if (std::abs(xl - r.x()) < SURFACE_COINCIDENT) {
    if (u.x() < 0.) {
      i--;
      on[0] = 1;
    } else {
      on[0] = -1;
    }
  } else if (std::abs(xh - r.x()) < SURFACE_COINCIDENT) {
    if (u.x() < 0.) {
      on[0] = 1;
    } else {
      i++;
      on[0] = -1;
    }
  }

  if (std::abs(yl - r.y()) < SURFACE_COINCIDENT) {
    if (u.y() < 0.) {
      j--;
      on[1] = 1;
    } else {
      on[1] = -1;
    }
  } else if (std::abs(yh - r.y()) < SURFACE_COINCIDENT) {
    if (u.y() < 0.) {
      on[1] = 1;
    } else {
      j++;
      on[1] = -1;
    }
  }

  if (std::abs(zl - r.z()) < SURFACE_COINCIDENT) {
    if (u.z() < 0.) {
      k--;
      on[2] = 1;
    } else {
      on[2] = -1;
    }
  } else if (std::abs(zh - r.z()) < SURFACE_COINCIDENT) {
    if (u.z() < 0.) {
      on[2] = 1;
    } else {
      k++;
      on[2] = -1;
    }
  }
}

std::pair<double, int> TrackLengthMeshTally::distance_to_next_index(
    const Position &r, const Direction &u, int i, int j, int k,
    const std::array<int, 3> &on) {
  // Get position at center of current tile
  double xc = r_low.x() + i * dx + 0.5 * dx;
  double yc = r_low.y() + j * dy + 0.5 * dy;
  double zc = r_low.z() + k * dz + 0.5 * dz;

  // Get relative position in cell
  Position r_tile(r.x() - xc, r.y() - yc, r.z() - zc);

  // Set our initial value for the distance and the index change
  double dist = INF;
  int key = 0;

  // Check all six sides
  double diff_xl = -dx * 0.5 - r_tile.x();
  double diff_xh = dx * 0.5 - r_tile.x();
  double diff_yl = -dy * 0.5 - r_tile.y();
  double diff_yh = dy * 0.5 - r_tile.y();
  double diff_zl = -dz * 0.5 - r_tile.z();
  double diff_zh = dz * 0.5 - r_tile.z();

  double d_xl = diff_xl / u.x();
  double d_xh = diff_xh / u.x();
  double d_yl = diff_yl / u.y();
  double d_yh = diff_yh / u.y();
  double d_zl = diff_zl / u.z();
  double d_zh = diff_zh / u.z();

  if (d_xl > 0. && d_xl < dist && on[0] != -1) {
    dist = d_xl;
    key = -1;
  }

  if (d_xh > 0. && d_xh < dist && on[0] != 1) {
    dist = d_xh;
    key = 1;
  }

  if (d_yl > 0. && d_yl < dist && on[1] != -1) {
    dist = d_yl;
    key = -2;
  }

  if (d_yh > 0. && d_yh < dist && on[1] != 1) {
    dist = d_yh;
    key = 2;
  }

  if (d_zl > 0. && d_zl < dist && on[2] != -1) {
    dist = d_zl;
    key = -3;
  }

  if (d_zh > 0. && d_zh < dist && on[2] != 1) {
    dist = d_zh;
    key = 3;
  }

  return {dist, key};
}

void TrackLengthMeshTally::update_indices(int key, int &i, int &j, int &k,
                                          std::array<int, 3> &on) {
  // Must initially fill with zero, so that we don't stay on top
  // of other surfaces the entire time
  on.fill(0);

  switch (key) {
    case -1:
      i--;
      on[0] = 1;
      break;

    case 1:
      i++;
      on[0] = -1;
      break;

    case -2:
      j--;
      on[1] = 1;
      break;

    case 2:
      j++;
      on[1] = -1;
      break;

    case -3:
      k--;
      on[2] = 1;
      break;

    case 3:
      k++;
      on[2] = -1;
      break;

    default:
      break;
  }
}
