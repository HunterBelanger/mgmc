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
#include <materials/const_material.hpp>
#include <plotting/slice_plot.hpp>
#include <utils/error.hpp>

//=========================================================
// Material Holders
std::map<uint32_t, uint32_t> material_id_to_indx;
std::vector<std::shared_ptr<Material>> materials;
std::vector<double> majorant_xs;

ConstMaterial::ConstMaterial(std::vector<double> i_Et, std::vector<double> i_Ea,
                             std::vector<double> i_Ef, std::vector<double> i_nu,
                             std::vector<std::vector<double>> i_chi,
                             std::vector<std::vector<double>> i_Es, uint32_t id,
                             std::string name)
    : Material{id, name},
      Et_(i_Et),
      Ea_(i_Ea),
      Ef_(i_Ef),
      nu_(i_nu),
      chi_(i_chi),
      Es_(i_Es),
      ngroups(),
      single_chi(true) {
  // Go through Ef, if not zero, make fissionable
  for (const double& el : Ef_) {
    if (el > 0.) {
      fissionable = true;
      break;
    }
  }

  ngroups = static_cast<int>(Et_.size());

  single_chi = true;
  if (chi_.size() != 1) single_chi = false;
}

double ConstMaterial::Et(Position, int E) const {
  if (E >= 0 and E < ngroups) {
    return Et_[E];
  } else {
    std::string mssg = std::to_string(E) + " is an invalide energy.";
    fatal_error(mssg, __FILE__, __LINE__);
    return 0.;
  }
}

double ConstMaterial::Ea(Position, int E) const {
  if (E >= 0 and E < ngroups) {
    return Ea_[E];
  } else {
    std::string mssg = std::to_string(E) + " is an invalide energy.";
    fatal_error(mssg, __FILE__, __LINE__);
    return 0.;
  }
}

double ConstMaterial::Ef(Position, int E) const {
  if (E >= 0 and E < ngroups) {
    return Ef_[E];
  } else {
    std::string mssg = std::to_string(E) + " is an invalide energy.";
    fatal_error(mssg, __FILE__, __LINE__);
    return 0.;
  }
}

double ConstMaterial::nu(Position, int E) const {
  if (E >= 0 and E < ngroups) {
    return nu_[E];
  } else {
    std::string mssg = std::to_string(E) + " is an invalide energy.";
    fatal_error(mssg, __FILE__, __LINE__);
    return 0.;
  }
}

const std::vector<double>& ConstMaterial::chi(Position, int E) const {
  if (single_chi)
    return chi_[0];
  else {
    if (E < 0 || E >= ngroups) {
      std::string mssg = std::to_string(E) + " is an invalide energy.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    return chi_[E];
  }
}

const std::vector<double>& ConstMaterial::Es(Position, int E) const {
  if (E < 0 || E >= ngroups) {
    std::string mssg = std::to_string(E) + " is an invalide energy.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  return Es_[E];
}

bool ConstMaterial::verify() const {
  // Check all dimensions
  if (ngroups != static_cast<int>(Ef_.size()) or
      ngroups != static_cast<int>(Ea_.size()) or
      ngroups != static_cast<int>(Es_.size()) or
      ngroups != static_cast<int>(nu_.size()) or
      ngroups != static_cast<int>(chi_.size()))
    return false;

  for (int i = 0; i < ngroups; i++) {
    if (ngroups != static_cast<int>(Es_[i].size())) return false;
  }

  // Check values are consistent
  for (int g = 0; g < ngroups; g++) {
    if (Et_[g] < 0 or Ea_[g] < 0 or Ef_[g] < 0 or nu_[g] < 0) return false;
    if (Ef_[g] > Ea_[g]) return false;

    // Get total scatter
    double Es_g = 0.0;
    for (int j = 0; j < ngroups; j++) {
      if (Es_[g][j] < 0.0)
        return false;
      else
        Es_g += Es_[g][j];
    }

    double err = std::abs(Et_[g] - (Ea_[g] + Es_g)) / Et_[g];
    if (err > 0.01) return false;
  }

  // Check chi
  if (chi_.size() != 1 && static_cast<int>(chi_.size()) != ngroups)
    return false;

  if (chi_.size() == 1) {
    if (static_cast<int>(chi_[0].size()) != ngroups)
      return false;
    else {
      for (const auto& val : chi_[0]) {
        if (val < 0.) return false;
      }
    }
  } else if (static_cast<int>(chi_.size()) == ngroups) {
    for (int i = 0; i < ngroups; i++) {
      if (static_cast<int>(chi_[i].size()) != ngroups)
        return false;
      else {
        for (int j = 0; j < ngroups; j++) {
          if (chi_[i][j] < 0.) return false;
        }
      }
    }
  }

  return true;
}

//=======================================================================
// Non-Member Functions
void make_const_material(YAML::Node mat_node) {
  // Get material id
  uint32_t id;
  if (mat_node["id"] && mat_node["id"].IsScalar()) {
    id = mat_node["id"].as<uint32_t>();
  } else {
    std::string mssg = "Material is missing a valid id.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get material name
  std::string name = "";
  if (mat_node["name"] && mat_node["name"].IsScalar()) {
    name = mat_node["name"].as<std::string>();
  }

  // Make vector for total XS
  std::vector<double> sig_tot;
  if (mat_node["total"] && mat_node["total"].IsSequence()) {
    for (size_t i = 0; i < mat_node["total"].size(); i++) {
      sig_tot.push_back(mat_node["total"][i].as<double>());
    }

    if (majorant_xs.size() == 0)
      majorant_xs = sig_tot;
    else {
      // Ensure proper number of elements
      if (sig_tot.size() != majorant_xs.size()) {
        std::string mssg = "Materials have different numbers of groups.";
        fatal_error(mssg, __FILE__, __LINE__);
      }

      // Get majorant XS
      for (size_t g = 0; g < majorant_xs.size(); g++) {
        if (sig_tot[g] > majorant_xs[g]) majorant_xs[g] = sig_tot[g];
      }
    }

  } else {
    std::string mssg = "Material is missing a valid total xs definition.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make vector for abs XS
  std::vector<double> sig_abs;
  if (mat_node["absorption"] && mat_node["absorption"].IsSequence()) {
    for (size_t i = 0; i < mat_node["absorption"].size(); i++) {
      sig_abs.push_back(mat_node["absorption"][i].as<double>());
    }

  } else {
    std::string mssg = "Material is missing a valid absorption xs definition.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make vector for fission XS
  std::vector<double> sig_fis;
  if (mat_node["fission"] && mat_node["fission"].IsSequence()) {
    for (size_t i = 0; i < mat_node["fission"].size(); i++) {
      sig_fis.push_back(mat_node["fission"][i].as<double>());
    }

  } else {
    std::string mssg = "Material is missing a valid fission xs definition.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make vector for nu
  std::vector<double> nu;
  if (mat_node["nu"] && mat_node["nu"].IsSequence()) {
    for (size_t i = 0; i < mat_node["nu"].size(); i++) {
      nu.push_back(mat_node["nu"][i].as<double>());
    }

  } else {
    std::string mssg = "Material is missing a valid nu definition.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make vector for chi
  /*std::vector<double> chi;
  if(mat_node["chi"] && mat_node["chi"].IsSequence()) {
    for(size_t i = 0; i < mat_node["chi"].size(); i++) {
      chi.push_back(mat_node["chi"][i].as<double>());
    }

  } else {
    std::string mssg = "Material is missing a valid fission xs definition.";
    fatal_error(mssg,__FILE__,__LINE__);
  }*/

  // Make chi matrix
  std::vector<std::vector<double>> chi;
  if (mat_node["chi"] && mat_node["chi"].IsSequence()) {
    for (size_t i = 0; i < mat_node["chi"].size(); i++) {
      // Make temp vector
      std::vector<double> tmp;
      if (mat_node["chi"][i].IsSequence()) {
        for (size_t j = 0; j < mat_node["chi"][i].size(); j++) {
          tmp.push_back(mat_node["chi"][i][j].as<double>());
        }
        chi.push_back(tmp);
        tmp.clear();
      } else {
        std::string mssg = "Material is missing valid chi matrix.";
        fatal_error(mssg, __FILE__, __LINE__);
      }
    }
  } else {
    std::string mssg = "Material is missing a valid scatter matrix.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make scattering matrix
  std::vector<std::vector<double>> sig_sct;
  if (mat_node["scatter"] && mat_node["scatter"].IsSequence()) {
    for (size_t i = 0; i < mat_node["scatter"].size(); i++) {
      // Make temp vector
      std::vector<double> tmp;
      if (mat_node["scatter"][i].IsSequence()) {
        for (size_t j = 0; j < mat_node["scatter"][i].size(); j++) {
          tmp.push_back(mat_node["scatter"][i][j].as<double>());
        }
        sig_sct.push_back(tmp);
        tmp.clear();
      } else {
        std::string mssg = "Material scatter matrix must be list of lists.";
        fatal_error(mssg, __FILE__, __LINE__);
      }
    }
  } else {
    std::string mssg = "Material is missing a valid scatter matrix.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Make sure material ID not yet taken
  if (material_id_to_indx.find(id) != material_id_to_indx.end()) {
    std::string mssg =
        "Material with id " + std::to_string(id) + "already exists.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  material_id_to_indx[id] = materials.size();

  // Get the color for the material if one is provided. Only needed for
  // plotting purposes
  plotter::Pixel mat_color;
  if (mat_node["color"] && mat_node["color"].IsSequence() &&
      mat_node["color"].size() == 3) {
    uint8_t R = mat_node["color"][0].as<int>();
    uint8_t G = mat_node["color"][1].as<int>();
    uint8_t B = mat_node["color"][2].as<int>();
    mat_color = plotter::Pixel(R, G, B);
    plotter::material_id_to_color[id] = mat_color;
  }

  materials.push_back(std::make_shared<ConstMaterial>(
      sig_tot, sig_abs, sig_fis, nu, chi, sig_sct, id, name));
}
