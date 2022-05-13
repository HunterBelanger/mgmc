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
#ifndef RECT_LATTICE_H
#define RECT_LATTICE_H

#include <geometry/lattice.hpp>

class RectLattice : public Lattice {
 public:
  RectLattice(uint32_t nx, uint32_t ny, uint32_t nz, double px, double py,
              double pz, double xl, double yl, double zl, uint32_t i_id,
              std::string i_name);
  ~RectLattice() = default;

  bool is_inside(Position r, Direction u) const override;

  std::array<int32_t, 3> get_tile(Position r, Direction u) const override;

  std::shared_ptr<Cell> get_cell(Position r, Direction u,
                                 int32_t on_surf) const override;

  Cell* get_cell_naked_ptr(Position r, Direction u,
                           int32_t on_surf) const override;

  std::shared_ptr<Cell> get_cell(std::vector<GeoLilyPad>& stack, Position r,
                                 Direction u, int32_t on_surf) const override;

  double distance_to_tile_boundary(Position r_local, Direction u,
                                   std::array<int32_t, 3> tile) const override;

  void set_elements(std::vector<int32_t> univs) override;

 private:
  // Values required to describe a rectilinear lattice
  uint32_t Nx, Ny, Nz;  // Number of elements along each axis
  double Px, Py, Pz;    // Pitch along each axis
  double Xl, Yl, Zl;    // Lowest value of lattice along each coordinate

  size_t linear_index(uint32_t nx, uint32_t ny, uint32_t nz) const;

  Position tile_center(int nx, int ny, int nz) const;

};  // RectLattice

void make_rect_lattice(YAML::Node latt_node, YAML::Node input);

#endif
