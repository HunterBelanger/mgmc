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
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <geometry/cell.hpp>
#include <geometry/geo_lily_pad.hpp>
#include <geometry/lattice.hpp>
#include <geometry/surfaces/surface.hpp>
#include <geometry/universe.hpp>
#include <simulation/cancel_bin.hpp>
#include <unordered_map>

struct CancelBinKey {
  size_t i, j, k;

  bool operator==(const CancelBinKey& other) const {
    return ((i == other.i) && (j == other.j) && (k == other.k));
  }
};

namespace geometry {

//==========================================================================
// Vectors for Geometry Objects
//--------------------------------------------------------------------------
// All surfaces in problem
extern std::vector<std::shared_ptr<Surface>> surfaces;
// Indicies to all surfaces which are vacuum or reflective
extern std::vector<uint32_t> boundaries;

// All cells in problem
extern std::vector<std::shared_ptr<Cell>> cells;

// Root Universe containing problem geometry, may be a
// CellUniverse or a LatticeUniverse
extern std::shared_ptr<Universe> root_universe;

// All universes in problem(including root universe at 0)
extern std::vector<std::shared_ptr<Universe>> universes;

// All lattices in problem
extern std::vector<std::shared_ptr<Lattice>> lattices;

// All cancelation bins for carter tracking
extern std::unordered_map<CancelBinKey, CancelBin> cancel_bins;
extern Position cancel_bins_low;
extern std::array<size_t, 3> cancel_bins_shape;
extern std::array<double, 3> cancel_bins_pitch;

//==========================================================================
// Boundary struct to contain information about boundary crossings
struct Boundary {
  Boundary(double d, std::shared_ptr<Surface> surf, BoundaryType bound);

  double distance;
  std::shared_ptr<Surface> surface;
  BoundaryType boundary_type;
  int32_t token;
};

//==========================================================================
// Functions
std::shared_ptr<Cell> get_cell(const Position& r, const Direction& u,
                               int32_t on_surf = 0);

std::shared_ptr<Cell> get_cell(std::vector<GeoLilyPad>& stack,
                               const Position& r, const Direction& u,
                               int32_t on_surf = 0);

Boundary get_boundary(const Position& r, const Direction& u,
                      int32_t on_surf = 0);

int32_t id_to_token(int32_t id);

void do_reflection(Particle& p, Boundary surface);

}  // namespace geometry

namespace std {
template <>
struct hash<CancelBinKey> {
  std::size_t operator()(const CancelBinKey& key) const {
    size_t Nz = geometry::cancel_bins_shape[2];
    size_t Ny = geometry::cancel_bins_shape[1];
    return key.k + Nz * (key.j + Ny * key.i);
  }
};
}  // namespace std

#endif