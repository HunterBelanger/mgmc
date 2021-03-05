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
#ifndef MG_CONST_MATERIAL_H
#define MG_CONST_MATERIAL_H

#include <yaml-cpp/yaml.h>

#include <map>
#include <materials/material.hpp>

// ConstMaterial class
class ConstMaterial : public Material {
 public:
  ConstMaterial(std::vector<double> i_Et, std::vector<double> i_Ea,
                std::vector<double> i_Ef, std::vector<double> i_nu,
                std::vector<std::vector<double>> i_chi,
                std::vector<std::vector<double>> i_Es, uint32_t id,
                std::string name = "");
  ~ConstMaterial() = default;

  double Et(Position r, int E) const;
  double Ea(Position r, int E) const;
  double Ef(Position r, int E) const;
  double nu(Position r, int E) const;
  const std::vector<double>& Es(Position r, int E) const;
  const std::vector<double>& chi(Position r, int E) const;
  bool verify() const;

  int num_groups() const { return ngroups; }

 private:
  std::vector<double> Et_;
  std::vector<double> Ea_;
  std::vector<double> Ef_;
  std::vector<double> nu_;
  std::vector<std::vector<double>> chi_;
  std::vector<std::vector<double>> Es_;
  int ngroups;
  bool single_chi;

};  // ConstMaterail

//===========================================================================
// Non-Member Functions
void make_const_material(YAML::Node mat_node);

#endif  // MG_CONST_MATERIAL_H