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
#ifndef MG_CANCEL_BIN_H
#define MG_CANCEL_BIN_H

#include <materials/material.hpp>
#include <simulation/particle.hpp>
#include <utils/rng.hpp>
#include <vector>

class CancelBin {
 public:
  CancelBin();
  CancelBin(std::shared_ptr<Material> mat);

  // Copy constructor does not actually do any copying. It simple
  // initializes a new, empty CancelBin. This is only here to
  // allow resize of vector of cancel bins, so we can re-shape the
  // global NDArray of CancelBins, and have a mutex to load
  // particles in parallel.
  CancelBin(const CancelBin&);

  CancelBin& operator=(const CancelBin& other);

  // Method to attempt to add particle to bin. Returns true if
  // particle was sucessfully added (meaning the particle was
  // located inside the bin) or false otherwise (meaning the
  // particle was outisde the bin).
  void add_particle(Particle& p, double Xl, double Xh, double Yl, double Yh,
                    double Zl, double Zh);

  void perform_cancelation(double Xl, double Xh, double Yl, double Yh,
                           double Zl, double Zh);

  std::vector<Particle> get_uniform_particles(std::shared_ptr<RNG> rng,
                                              double Xl, double Xh, double Yl,
                                              double Yh, double Zl, double Zh);

 private:
  std::vector<Particle*> bin;
  double uniform_weight = 0.;
  std::shared_ptr<Material> material;

  bool is_inside(const Position& r, double Xl, double Xh, double Yl, double Yh,
                 double Zl, double Zh);
  double get_beta(Position r_parent, double Esmp, double Ef, double Xl,
                  double Xh, double Yl, double Yh, double Zl, double Zh);
  double get_f(Position r, Position r_parent, double Esmp, double Ef);

};  // CancelBin

#endif  // MG_CANCEL_BIN_H
