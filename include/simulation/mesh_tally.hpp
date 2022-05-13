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
#ifndef MESH_TALLY_H
#define MESH_TALLY_H

#include <yaml-cpp/yaml.h>

#include <materials/material_helper.hpp>
#include <materials/nuclide.hpp>
#include <memory>
#include <ndarray.hpp>
#include <simulation/particle.hpp>

class MeshTally {
 public:
  enum class TallyQuantity {
    Flux,
    Total,
    Elastic,
    Absorption,
    Fission,
    MT,
    RealFlux,
    ImgFlux,
    MagFlux,
    MagSqrFlux,
    RealTotal,
    ImgTotal,
    RealElastic,
    ImgElastic,
    RealAbsorption,
    ImgAbsorption,
    RealFission,
    ImgFission,
    RealMT,
    ImgMT
  };

  MeshTally(Position low, Position hi, uint64_t nx, uint64_t ny, uint64_t nz,
            const std::vector<double> &ebounds, std::string fname);
  virtual ~MeshTally() = default;

  void set_net_weight(double W);

  // multiplier must be the same across ALL MPI processes, as it is
  // applied after the MPI reduction of the generation scores.
  void record_generation(double multiplier = 1.);

  void clear_generation();

  void write_tally();

 protected:
  Position r_low, r_hi;
  uint64_t Nx, Ny, Nz, g;
  double dx, dy, dz, net_weight;
  std::vector<double> energy_bounds;
  std::string fname;

  NDArray<double> tally_gen;
  NDArray<double> tally_avg;
  NDArray<double> tally_var;

  virtual std::string quantity_str() const = 0;
};

#endif
