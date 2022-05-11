/*    Author(s): Vincent Rouvreau
 *
 *    Copyright (C) 2021 Inria - MIT license
 *
 *    Modification(s):
 */

#include <iostream>
#include <cmath>  // for cos, sin
#include <string>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Implicit_manifold_intersection_oracle.h>  // for Gudhi::coxeter_triangulation::make_oracle
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Cell_complex.h>

#include <gudhi/IO/build_mesh_from_cell_complex.h>
#include <gudhi/IO/output_meshes_to_medit.h>

#include <Eigen/Dense>

using namespace Gudhi::coxeter_triangulation;

const double pi = std::acos(-1);

struct Function_knot {
  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(k_);

    double x = p(0);
    double x_2 = x*x;
    double y = p(1);
    double y_2 = y*y;
    double z = p(2);
    double z_2 = z*z;
    double t = p(3);
    double t_2 = t*t;
    double rho = x_2 + y_2 + (z_2 - t_2) * (z_2 - t_2) + 4. * z_2 * t_2;

    result(0) = x_2 + y_2 + z_2 + t_2 - 1.;
    double f_1 = (z_2 - t_2) * rho + x * (8. * x_2 - 2. * rho);
    double f_2 = 2. * std::sqrt(2) * x * z * t + y * (8. * x_2 - rho);
    result(1) = c_ * f_1 + s_*f_2;
    return result;
  }

  /** \brief Returns the domain dimension. Same as the ambient dimension of the knot. */
  std::size_t amb_d() const { return d_; };

  /** \brief Returns the codomain dimension. Same as the codimension of the knot. */
  std::size_t cod_d() const { return k_; };

  /** \brief Returns a point on the knot. */
  Eigen::VectorXd seed() const {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(d_);
    result(2) = 1./std::sqrt(2.);
    result(3) = 1./std::sqrt(2.);
    return result;
  }

  Function_knot(double cosinus, double sinus)
      : c_(cosinus), s_(sinus), m_(1), k_(2), d_(4) {}

 private:
  double c_, s_;
  std::size_t m_, k_, d_;
};

struct Function_fibration {
  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(k_);

    double x = p(0);
    double x_2 = x*x;
    double y = p(1);
    double y_2 = y*y;
    double z = p(2);
    double z_2 = z*z;
    double t = p(3);
    double t_2 = t*t;
    double rho = x_2 + y_2 + (z_2 - t_2) * (z_2 - t_2) + 4. * z_2 * t_2;

    double f_1 = (z_2 - t_2) * rho + x * (8. * x_2 - 2. * rho);
    double f_2 = 2. * std::sqrt(2) * x * z * t + y * (8. * x_2 - rho);
    // s*f1-c*f2>0, in other words:
    result(0) = c_ * f_2 - s_ * f_1;
    return result;
  }

  /** \brief Returns the domain dimension. Same as the ambient dimension of the knot. */
  std::size_t amb_d() const { return d_; };

  /** \brief Returns the codomain dimension. Same as the codimension of the knot. */
  std::size_t cod_d() const { return k_; };

  /** \brief Returns a point on the knot. */
  Eigen::VectorXd seed() const {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(d_);
    result(2) = 1./std::sqrt(2.);
    result(3) = 1./std::sqrt(2.);
    return result;
  }

  Function_fibration(double cosinus, double sinus)
      : c_(cosinus), s_(sinus), m_(1), k_(1), d_(4) {}

 private:
  double c_, s_;
  std::size_t m_, k_, d_;
};


int main(int argc, char** argv) {
  for (double angle = 0.; angle < 360.; angle +=5.) {
    double alpha = 2 * pi * (angle / 360.);
    std::cout << std::cos(alpha) << ", " << std::sin(alpha) << std::endl;
    // Oracle is a circle of radius 1
    auto oracle = make_oracle(Function_knot(std::cos(alpha), std::sin(alpha)),
                              Function_fibration(std::cos(alpha), std::sin(alpha))); 

    // Define a Coxeter triangulation.
    Coxeter_triangulation<> cox_tr(oracle.amb_d());
    // Theory forbids that a vertex of the triangulation lies exactly on the circle.
    // Add some offset to avoid algorithm degeneracies.
    cox_tr.change_offset(-Eigen::VectorXd::Random(oracle.amb_d()));
    // For a better manifold approximation, one can change the circle radius value or change the linear transformation
    // matrix.
    // The number of points and edges will increase with a better resolution.
    cox_tr.change_matrix(0.05 * cox_tr.matrix());

    // Manifold tracing algorithm
    using Out_simplex_map = typename Manifold_tracing<Coxeter_triangulation<> >::Out_simplex_map;

    std::vector<Eigen::VectorXd> seed_points(1, oracle.seed());
    Out_simplex_map interior_simplex_map, boundary_simplex_map;
    manifold_tracing_algorithm(seed_points, cox_tr, oracle, interior_simplex_map, boundary_simplex_map);

    // Constructing the cell complex
    std::size_t intr_d = oracle.amb_d() - oracle.cod_d();
    Cell_complex<Out_simplex_map> cell_complex(intr_d);
    cell_complex.construct_complex(interior_simplex_map, boundary_simplex_map);

    // List of Hasse_cell pointers to retrieve vertices values from edges
    std::map<Cell_complex<Out_simplex_map>::Hasse_cell*, std::size_t> vi_map;

    auto mesh = build_mesh_from_cell_complex(cell_complex,
                                             Configuration(true, true, true, 1, 5, 3),
                                             Configuration(true, true, true, 2, 13, 14));
  
    // Output the cell complex to a file readable by medit
    std::string meshfilename("seifert_surface_");
    meshfilename += std::to_string(angle);
    output_meshes_to_medit(4, meshfilename, mesh);
  }
  return 0;
}
