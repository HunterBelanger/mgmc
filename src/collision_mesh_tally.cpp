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
#include <sstream>
#include <utils/error.hpp>

void CollisionMeshTally::score_collision(const Particle &p,
                                         MaterialHelper &mat) {
  double Et = mat.Et(p.E());

  double scr = 1. / (Et * net_weight);

  /*if (p.wgt() == 0. && p.wgt2() == 0.) {
    std::stringstream out;
    out << " Zero collision weights.\n";
    std::cout << out.str();
  }*/

  int i = std::floor((p.r().x() - r_low.x()) / dx);
  int j = std::floor((p.r().y() - r_low.y()) / dy);
  int k = std::floor((p.r().z() - r_low.z()) / dz);
  int l = -1;

  // Get energy index with linear search
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

  if (i >= 0 && i < static_cast<int>(Nx) && j >= 0 &&
      j < static_cast<int>(Ny) && k >= 0 && k < static_cast<int>(Nz)) {
    uint64_t ui = static_cast<uint64_t>(i);
    uint64_t uj = static_cast<uint64_t>(j);
    uint64_t uk = static_cast<uint64_t>(k);
    uint64_t uE = static_cast<uint64_t>(l);

    // Must multiply score by correct factor depending on the observed
    // quantity. q = 0 corresponds to just the flux, so no modification.
    switch (quantity) {
      case TallyQuantity::Flux:
        scr *= p.wgt();
        break;

      case TallyQuantity::Elastic:
        scr *= p.wgt() * mat.Eelastic(p.E());
        break;

      case TallyQuantity::Absorption:
        scr *= p.wgt() * mat.Ea(p.E());
        break;

      case TallyQuantity::Fission:
        scr *= p.wgt() * mat.Ef(p.E());
        break;

      case TallyQuantity::Total:
        scr *= p.wgt() * Et;
        break;

      case TallyQuantity::MT:
        scr *= p.wgt() * mat.Emt(mt, p.E());
        break;

      case TallyQuantity::RealFlux:
        scr *= p.wgt();
        break;

      case TallyQuantity::ImgFlux:
        scr *= p.wgt2();
        break;

      case TallyQuantity::RealTotal:
        scr *= p.wgt() * Et;
        break;

      case TallyQuantity::ImgTotal:
        scr *= p.wgt2() * Et;
        break;

      case TallyQuantity::RealElastic:
        scr *= p.wgt() * mat.Eelastic(p.E());
        break;

      case TallyQuantity::ImgElastic:
        scr *= p.wgt2() * mat.Eelastic(p.E());
        break;

      case TallyQuantity::RealAbsorption:
        scr *= p.wgt() * mat.Ea(p.E());
        break;

      case TallyQuantity::ImgAbsorption:
        scr *= p.wgt2() * mat.Ea(p.E());
        break;

      case TallyQuantity::RealFission:
        scr *= p.wgt() * mat.Ef(p.E());
        break;

      case TallyQuantity::ImgFission:
        scr *= p.wgt2() * mat.Ef(p.E());
        break;

      case TallyQuantity::RealMT:
        scr *= p.wgt() * mat.Emt(mt, p.E());
        break;

      case TallyQuantity::ImgMT:
        scr *= p.wgt2() * mat.Emt(mt, p.E());
        break;

      case TallyQuantity::MagFlux:
        scr *= std::sqrt(p.wgt() * p.wgt() + p.wgt2() * p.wgt2());
        break;

      case TallyQuantity::MagSqrFlux:
        scr *= (p.wgt() * p.wgt() + p.wgt2() * p.wgt2());
        break;
    }
#ifdef _OPENMP
#pragma omp atomic
#endif
    tally_gen(uE, ui, uj, uk) += scr;
  }
}
