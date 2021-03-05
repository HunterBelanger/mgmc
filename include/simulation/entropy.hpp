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
#ifndef MG_ENTROPY_H
#define MG_ENTROPY_H

#include <array>
#include <cmath>
#include <utils/position.hpp>
#include <vector>

class Entropy {
 public:
  enum Sign { Positive, Negative, Total };

 public:
  Entropy(Position low, Position up, std::array<uint32_t, 3> shp, Sign sgn)
      : lower_corner{low},
        upper_corner{up},
        shape{shp},
        nbins{0},
        dx{},
        dy{},
        dz{},
        bins{},
        total_weight{},
        sign{sgn} {
    nbins = shape[0] * shape[1] * shape[2];
    bins.resize(nbins, 0.);
    total_weight = 0.;
    dx = (upper_corner.x() - lower_corner.x()) / static_cast<double>(shape[0]);
    dy = (upper_corner.y() - lower_corner.y()) / static_cast<double>(shape[1]);
    dz = (upper_corner.z() - lower_corner.z()) / static_cast<double>(shape[2]);
  }
  ~Entropy() = default;

  void add_point(const Position& r, const double& w) {
    // Get bin idecies
    int32_t nx =
        static_cast<int32_t>(std::floor((r.x() - lower_corner.x()) / dx));
    int32_t ny =
        static_cast<int32_t>(std::floor((r.y() - lower_corner.y()) / dy));
    int32_t nz =
        static_cast<int32_t>(std::floor((r.z() - lower_corner.z()) / dz));

    // Make sure position is in a valid bin
    if (nx >= 0 && nx < static_cast<int32_t>(shape[0]) && ny >= 0 &&
        ny < static_cast<int32_t>(shape[1]) && nz >= 0 &&
        nz < static_cast<int32_t>(shape[2])) {
// Add to total weight, no matter sign being tallied
#pragma omp atomic
      this->total_weight += w;

      // Get sign of particle
      Sign p_sign;
      if (w > 0.)
        p_sign = Sign::Positive;
      else
        p_sign = Sign::Negative;

      if (sign == Sign::Total) p_sign = Sign::Total;

      // Add weight to bin if sign agrees
      if (p_sign == sign) {
#pragma omp atomic
        this->bins[(shape[1] * shape[2]) * nx + (shape[2]) * ny + nz] += w;
      }

    } else
      std::cout << " Missing entropy particle at " << r << "\n";
  }

  double calculate_entropy() const {
    // Set sum to zero
    double sum = 0.;

    /*for(const auto& b : bins) {
      double p = std::fabs(b) / total_weight;
      if(p != 0.)
        sum -= p*std::log2(p);
    }*/

#pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < this->bins.size(); i++) {
      double p = std::fabs(this->bins[i]) / this->total_weight;

      if (p > 1.0) {
#pragma omp critical
        {
          std::cout << " Negative entropy: p = " << p
                    << ", bin = " << this->bins[i]
                    << ", total_weight = " << this->total_weight << "\n";
        }
      }

      if (p != 0.) {
#pragma omp atomic
        sum -= p * std::log2(p);
      }
    }

    return sum;
  }

  void zero() {
    // Zero all bins
    for (auto& b : this->bins) b = 0.;

    this->total_weight = 0.;
  }

  double total_wgt() const { return total_weight; }

 private:
  Position lower_corner;
  Position upper_corner;
  std::array<uint32_t, 3> shape;
  uint32_t nbins;
  double dx, dy, dz;
  std::vector<double> bins;
  double total_weight;
  Sign sign;

};  // Etnropy

#endif  // MG_ENTROPY_H
