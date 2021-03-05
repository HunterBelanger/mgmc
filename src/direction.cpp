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
#include <utils/direction.hpp>

Direction::Direction(double i_x, double i_y, double i_z)
    : Vector{i_x, i_y, i_z} {
  double magnitude = this->norm();
  this->x_ /= magnitude;
  this->y_ /= magnitude;
  this->z_ /= magnitude;
}

//============================================================================
// Addition Operators
Vector operator+(const Direction& d1, const Direction& d2) {
  return Vector(d1.x() + d2.x(), d1.y() + d2.y(), d1.z() + d2.z());
}

Vector operator+(const Direction& d, const Vector& v) {
  return Vector(d.x() + v.x(), d.y() + v.y(), d.z() + v.z());
}

Vector operator+(const Vector& v, const Direction& d) {
  return Vector(d.x() + v.x(), d.y() + v.y(), d.z() + v.z());
}

//============================================================================
// Subtraction Operators
Vector operator-(const Direction& d1, const Direction& d2) {
  return Vector(d1.x() - d2.x(), d1.y() - d2.y(), d1.z() - d2.z());
}

Vector operator-(const Direction& d, const Vector& v) {
  return Vector(d.x() - v.x(), d.y() - v.y(), d.z() - v.z());
}

Vector operator-(const Vector& v, const Direction& d) {
  return Vector(v.x() - d.x(), v.y() - d.y(), v.z() - d.z());
}

//============================================================================
// Dot Product Operators
double operator*(const Direction& d, const Direction& v) { return d.dot(v); }

double operator*(const Direction& d, const Vector& v) { return d.dot(v); }

double operator*(const Vector& v, const Direction& d) { return d.dot(v); }

//============================================================================
// Scaling Operators
Vector operator*(const Direction& d, double c) {
  return Vector(d.x() * c, d.y() * c, d.z() * c);
}

Vector operator*(double c, const Direction& d) {
  return Vector(d.x() * c, d.y() * c, d.z() * c);
}

Vector operator/(const Direction& d, double c) {
  return Vector(d.x() / c, d.y() / c, d.z() / c);
}

std::ostream& operator<<(std::ostream& output, const Direction& d) {
  output << "<<" << d.x() << "," << d.y() << "," << d.z() << ">>";
  return output;
}