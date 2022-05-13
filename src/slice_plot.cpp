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
#include <fstream>
#include <geometry/cell.hpp>
#include <geometry/geometry.hpp>
#include <memory>
#include <plotting/slice_plot.hpp>
#include <utils/error.hpp>

namespace plotter {

SlicePlot::SlicePlot(std::string fname, uint64_t pwidth, uint64_t pheight,
                     double width, double height, Position origin, Basis basis,
                     ColorBy colorby, Pixel background)
    : file_name_{fname},
      plot_width_{pwidth},
      plot_height_{pheight},
      basis_{basis},
      colorby_{colorby},
      origin_{origin},
      background_{background},
      image_matrix(),
      width_{width},
      height_{height},
      rng(2617257382),
      create_color_mutex{} {}

// Other Methods
void SlicePlot::generate_plot() {
  // Calculate number of elements
  uint64_t n_pixels = plot_width_ * plot_height_;

  // Shape image matrix accordingly
  image_matrix.resize(n_pixels);

// Go through each pixel
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (size_t i = 0; i < plot_height_; i++) {
    for (size_t j = 0; j < plot_width_; j++) {
      // Get cell pointer for pixel location
      Position pixel_r = get_pixel_position(i, j);
      std::shared_ptr<::Cell> cell =
          geometry::root_universe->get_cell(pixel_r, {1., 0., 0.}, 0);

      uint64_t indx = i * plot_width_ + j;
      Pixel pixel_color = background_;
      // If pointer isn't nullpntr, get pixel
      if (cell != nullptr) {
        if (colorby_ == ColorBy::Cell) {
          // Do same check twice with mutex to make thread safe
          if (cell_id_to_color.find(cell->id()) == cell_id_to_color.end()) {
            // Check if cell id is in id_to_pixel
            create_color_mutex.lock();
            if (cell_id_to_color.find(cell->id()) == cell_id_to_color.end()) {
              // Get new random color for id
              cell_id_to_color[cell->id()] = get_random_color();
            }
            create_color_mutex.unlock();
          }
          pixel_color = cell_id_to_color[cell->id()];
        } else {
          // TODO Color by material
          std::shared_ptr<::Material> material = cell->material();
          // Do same check twice with mutex to make thread safe
          if (material_id_to_color.find(material->id()) ==
              material_id_to_color.end()) {
            // Check if cell id is in id_to_pixel
            create_color_mutex.lock();
            if (material_id_to_color.find(material->id()) ==
                material_id_to_color.end()) {
              // Get new random color for id
              material_id_to_color[material->id()] = get_random_color();
            }
            create_color_mutex.unlock();
          }
          pixel_color = material_id_to_color[material->id()];
        }
      }
      // Set pixel
      image_matrix[indx] = pixel_color;
    }  // For j
  }    // For i which is parallel
}

void SlicePlot::write() const {
  if (image_matrix.size() > 0) {
    std::ofstream file(file_name_ + ".ppm");
    file << "P6\n";

    // Get width and height as string
    std::string width = std::to_string(plot_width_);
    std::string height = std::to_string(plot_height_);
    file << width << ' ' << height << '\n';

    // Set max value
    file << "255\n";

    // Write all pixels from matrix. First get pointer
    const char *img_data = reinterpret_cast<const char *>(image_matrix.data());
    // Then get number of bytes. Use sizeof to be sure works on all systems
    size_t nbytes = image_matrix.size() * sizeof(Pixel);
    // Write instead of streaming seems to be much faster
    file.write(img_data, nbytes);

    // Close file
    file.close();
  }
}

Position SlicePlot::get_pixel_position(uint64_t i, uint64_t j) const {
  // Make sure indicies are valid. i goes down so is height, j goes
  // across so is width
  if (i >= plot_height_ || j >= plot_width_) {
    std::string mssg = "Trying to deffine invalid pixel for plot.";
    fatal_error(mssg, __FILE__, __LINE__);
  }

  if (basis_ == Basis::XY) {
    // x is i going down, j is y going across
    // First get height and width of a pixel
    double dx = height_ / static_cast<double>(plot_height_);
    double dy = width_ / static_cast<double>(plot_width_);

    // Get upper corner off plot
    double x_low = origin_.x() - (0.5 * height_);
    double y_low = origin_.y() - (0.5 * width_);

    // Get coordinate of pixel
    double x = (static_cast<double>(i) + 0.5) * dx + x_low;
    double y = (static_cast<double>(j) + 0.5) * dy + y_low;

    // Return coordinate
    return {x, y, origin_.z()};
  } else if (basis_ == Basis::YZ) {
    // y is i going down, j is z going across
    // First get height and width of a pixel
    double dy = height_ / static_cast<double>(plot_height_);
    double dz = width_ / static_cast<double>(plot_width_);

    // Get upper corner off plot
    double y_low = origin_.y() - (0.5 * height_);
    double z_low = origin_.z() - (0.5 * width_);

    // Get coordinate of pixel
    double y = (static_cast<double>(i) + 0.5) * dy + y_low;
    double z = (static_cast<double>(j) + 0.5) * dz + z_low;

    // Return coordinate
    return {origin_.x(), y, z};
  } else {
    // x is i going down, j is z going across
    // First get height and width of a pixel
    double dx = height_ / static_cast<double>(plot_height_);
    double dz = width_ / static_cast<double>(plot_width_);

    // Get upper corner off plot
    double x_low = origin_.x() - (0.5 * height_);
    double z_low = origin_.z() - (0.5 * width_);

    // Get coordinate of pixel
    double x = (static_cast<double>(i) + 0.5) * dx + x_low;
    double z = (static_cast<double>(j) + 0.5) * dz + z_low;

    // Return coordinate
    return {x, origin_.y(), z};
  }
}

Pixel SlicePlot::get_random_color() {
  uint8_t r = static_cast<uint8_t>(255.0 * RNG::rand(rng));
  uint8_t g = static_cast<uint8_t>(255.0 * RNG::rand(rng));
  uint8_t b = static_cast<uint8_t>(255.0 * RNG::rand(rng));

  return {r, g, b};
}

}  // namespace plotter
