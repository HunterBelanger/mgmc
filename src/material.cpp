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
#include <materials/material.hpp>
#include <materials/mg_nuclide.hpp>
#include <memory>
#include <plotting/slice_plot.hpp>
#include <utils/error.hpp>
#include <utils/output.hpp>
#include <utils/settings.hpp>
#include <vector>

std::map<uint32_t, std::shared_ptr<Material>> materials;

void fill_mg_material(const YAML::Node &mat,
                      std::shared_ptr<Material> material) {
  std::string read_mssg = " Reading material " + material->get_name() + ".\n";
  Output::instance()->write(read_mssg);

  // Get total xs
  if (!mat["total"] || !mat["total"].IsSequence()) {
    std::string mssg = "Multigroup material missing total cross section.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  std::vector<double> Et = mat["total"].as<std::vector<double>>();

  // Get absorption xs
  if (!mat["absorption"] || !mat["absorption"].IsSequence()) {
    std::string mssg = "Multigroup material missing absorption cross section.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  std::vector<double> Ea = mat["absorption"].as<std::vector<double>>();

  // Get fission xs
  std::vector<double> Ef(settings::ngroups, 0.);
  if (mat["fission"] && mat["fission"].IsSequence()) {
    Ef.clear();
    Ef = mat["fission"].as<std::vector<double>>();
  } else if (mat["fission"]) {
    std::string mssg = "Multigroup material has invalid fission cross section.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get chi
  std::vector<double> chi(settings::ngroups, 0.);
  if (mat["chi"] && mat["chi"].IsSequence()) {
    chi.clear();
    chi = mat["chi"].as<std::vector<double>>();
  } else if (mat["chi"]) {
    std::string mssg = "Multigroup material has invalid chi.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get scatter xs
  if (!mat["scatter"] || !mat["scatter"].IsSequence()) {
    std::string mssg = "Multigroup material has invalid scattering matrix.";
    fatal_error(mssg, __FILE__, __LINE__);
  }
  std::vector<std::vector<double>> Es;

  for (size_t g = 0; g < mat["scatter"].size(); g++) {
    if (!mat["scatter"][g].IsSequence()) {
      std::string mssg = "Multigroup material has invalid scattering matrix.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    std::vector<double> Es_g = mat["scatter"][g].as<std::vector<double>>();
    Es.push_back(Es_g);
  }

  // Get nu
  std::vector<double> nu(settings::ngroups, 0.);
  std::vector<double> nu_delayed(settings::ngroups, 0.);
  if (mat["nu"]) {
    if (mat["nu"] && mat["nu"].IsSequence()) {
      nu.clear();
      nu = mat["nu"].as<std::vector<double>>();
    } else if (mat["nu"]) {
      std::string mssg = "Multigroup material has invalid nu.";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  } else if ((mat["nu_prompt"] && mat["nu_prompt"].IsSequence()) &&
             (mat["nu_delayed"] && mat["nu_delayed"].IsSequence())) {
    nu.clear();
    nu = mat["nu_prompt"].as<std::vector<double>>();

    nu_delayed.clear();
    nu_delayed = mat["nu_delayed"].as<std::vector<double>>();
  } else {
    std::string mssg = "Multigroup material missing valid nu entry.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get delayed group info
  std::vector<double> delayed_group_probs;
  std::vector<double> delayed_group_constants;
  if (mat["nu_prompt"] && mat["delayed_groups"] &&
      mat["delayed_groups"].IsMap()) {
    if (!mat["delayed_groups"]["probabilities"] ||
        !mat["delayed_groups"]["probabilities"].IsSequence()) {
      std::string mssg = "Material missing delayed group probabilities.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    if (!mat["delayed_groups"]["constants"] ||
        !mat["delayed_groups"]["constants"].IsSequence()) {
      std::string mssg = "Material missing delayed group constants.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    delayed_group_probs =
        mat["delayed_groups"]["probabilities"].as<std::vector<double>>();
    delayed_group_constants =
        mat["delayed_groups"]["constants"].as<std::vector<double>>();

    // Check all probs and constants now
    double prob_sum = 0.;
    for (const auto &v : delayed_group_probs) {
      if (v <= 0.) {
        std::string mssg = "Negative or zero delayed group probability found.";
        fatal_error(mssg, __FILE__, __LINE__);
      }

      prob_sum += v;
    }

    // Normalize probabilities
    for (auto &v : delayed_group_probs) v /= prob_sum;

    // Ensure proper number of constants
    if (delayed_group_probs.size() != delayed_group_constants.size()) {
      std::string mssg =
          "Different number of delayed group probabilities and group "
          "constants.";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Make sure no negative constants
    for (const auto &v : delayed_group_constants) {
      if (v <= 0.) {
        std::string mssg = "Negative or zero delayed group constant found.";
        fatal_error(mssg, __FILE__, __LINE__);
      }
    }
  } else if (!mat["nu"]) {
    std::string mssg = "Multigroup material missing delayed group info.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  //============================================================================
  // Verify XS data
  // Check all dimensions
  uint32_t ngroups = settings::ngroups;
  if (ngroups != Ef.size() || ngroups != Ea.size() || ngroups != Es.size() ||
      ngroups != chi.size() || ngroups != nu.size() ||
      (mat["nu_prompt"] && ngroups != nu_delayed.size())) {
    std::string mssg = "Multigroup material " + material->get_name() + " has ";
    mssg += "inconsisten number of entries.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  for (uint32_t i = 0; i < ngroups; i++) {
    if (ngroups != Es[i].size()) {
      std::string mssg =
          "Multigroup material " + material->get_name() + " has ";
      mssg += "non-square scatter matrix.";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  }

  // Check values are consistent
  for (uint32_t g = 0; g < ngroups; g++) {
    if (Et[g] < 0 || Ea[g] < 0 || Ef[g] < 0 || nu[g] < 0) {
      std::string mssg =
          "Multigroup material " + material->get_name() + " has ";
      mssg += "negative data.";
      fatal_error(mssg, __FILE__, __LINE__);
    }
    if (Ef[g] > Ea[g]) {
      std::string mssg = "Multigroup material " + material->get_name() + " Ef ";
      mssg += "is larger than Ea for group " + std::to_string(g) + ".";
      fatal_error(mssg, __FILE__, __LINE__);
    }

    // Get total scatter
    double Es_g = 0.0;
    for (size_t j = 0; j < ngroups; j++) {
      if (Es[g][j] < 0.0) {
        std::string mssg =
            "Multigroup material " + material->get_name() + " has ";
        mssg += "negative scatter matrix exlement.";
        fatal_error(mssg, __FILE__, __LINE__);
      } else {
        Es_g += Es[g][j];
      }
    }

    double err = std::abs(Et[g] - (Ea[g] + Es_g)) / Et[g];
    if (err > 0.001) {
      std::string mssg = "Multigroup material " + material->get_name() + " is ";
      mssg += "not self consistent.";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  }

  // Check chi
  if (chi.size() != ngroups) {
    std::string mssg = "Multigroup material " + material->get_name() + " has ";
    mssg += "chi of improper size.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  for (uint32_t j = 0; j < ngroups; j++) {
    if (chi[j] < 0.) {
      std::string mssg =
          "Multigroup material " + material->get_name() + " has ";
      mssg += "negative chi exlement.";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  }

  // Can add nuclide to material
  std::shared_ptr<Nuclide> nuclide;

  if (mat["nu"] && settings::mode != settings::SimulationMode::NOISE) {
    nuclide = std::make_shared<MGNuclide>(settings::energy_bounds, Et, Ea, Ef,
                                          nu, chi, Es);
  } else {
    nuclide = std::make_shared<MGNuclide>(
        settings::energy_bounds, Et, Ea, Ef, nu, nu_delayed, chi,
        delayed_group_probs, delayed_group_constants, Es);
  }

  material->add_component(1., nuclide);
}

void make_material(const YAML::Node &mat, const YAML::Node & /*xsdir*/,
                   bool plotting_mode) {
  std::shared_ptr<Material> material = std::make_shared<Material>();

  // Get material id
  if (mat["id"] && mat["id"].IsScalar()) {
    material->id_ = mat["id"].as<uint32_t>();
  } else {
    std::string mssg = "Material is missing a valid id.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // make sure material id isn't taken
  if (materials.find(material->id()) != materials.end()) {
    std::string mssg = "A material with id " + std::to_string(material->id()) +
                       " already exists.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  // Get material name
  if (mat["name"] && mat["name"].IsScalar()) {
    material->name_ = mat["name"].as<std::string>();
  }

  // Get the color for the material if one is provided. Only needed for
  // plotting purposes
  plotter::Pixel mat_color;
  if (mat["color"] && mat["color"].IsSequence() && mat["color"].size() == 3) {
    uint8_t R = mat["color"][0].as<int>();
    uint8_t G = mat["color"][1].as<int>();
    uint8_t B = mat["color"][2].as<int>();
    mat_color = plotter::Pixel(R, G, B);
    plotter::material_id_to_color[material->id()] = mat_color;
  }

  if (!plotting_mode) {
    if (settings::energy_mode == settings::EnergyMode::MG) {
      fill_mg_material(mat, material);
    } else if (settings::energy_mode == settings::EnergyMode::CE) {
      std::string mssg = "Continuous Energy Mode not supported.";
      fatal_error(mssg, __FILE__, __LINE__);
    } else {
      std::string mssg = "Unknown energy mode.";
      fatal_error(mssg, __FILE__, __LINE__);
    }
  }

  // Save material
  materials[material->id()] = material;
}
