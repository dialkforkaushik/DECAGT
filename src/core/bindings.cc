#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include "discrete_exterior_calculus.h"
#include "geometry.h"
#include "finite_element_exterior_calculus.h"
#include "definitions.h"
#include "simplicial_complex.h"
#include "core_utils.h"

namespace py = pybind11;


PYBIND11_MODULE(decagt, m) {
    py::class_<SimplicialComplex>(m, "SimplicialComplex")
        .def(py::init<>())
        .def(py::init<Vector2D &, Vector2I &>())
        .def("read_files", &SimplicialComplex::read_files)
        .def("compute_simplices", &SimplicialComplex::compute_simplices)
        .def("compute_boundary_matrices", &SimplicialComplex::compute_boundary_matrices, py::call_guard<py::gil_scoped_release>())
        .def("compute_adjacency1d", &SimplicialComplex::compute_adjacency1d)
        .def("compute_adjacency2d", &SimplicialComplex::compute_adjacency2d)
        .def("compute_elements", &SimplicialComplex::compute_elements)
        .def("compute_simplex_sub_simplices", &SimplicialComplex::compute_simplex_sub_simplices)
        .def("build_complex", &SimplicialComplex::build_complex)
        .def("circumcenter", &SimplicialComplex::circumcenter)
        .def_readonly("vertices", &SimplicialComplex::vertices)
        .def_readonly("simplex", &SimplicialComplex::simplex)
        .def_readonly("complex_dimension", &SimplicialComplex::complex_dimension)
        .def_readonly("embedding_dimension", &SimplicialComplex::embedding_dimension)
        .def_readonly("num_simplices", &SimplicialComplex::num_simplices)
        .def_readonly("boundary_matrices", &SimplicialComplex::boundary_matrices)
        .def_readonly("simplices", &SimplicialComplex::simplices)
        .def_readonly("elements", &SimplicialComplex::elements)
        .def_readonly("simplex_sub_simplices", &SimplicialComplex::simplex_sub_simplices)
        .def_readonly("adjacency1d", &SimplicialComplex::adjacency1d)
        .def_readonly("adjacency2d", &SimplicialComplex::adjacency2d);


    py::class_<GeometryComplex, SimplicialComplex>(m, "GeometryComplex")
        .def(py::init<>())
        .def(py::init<SimplicialComplex>())
        .def("compute_dual_volumes", &GeometryComplex::compute_dual_volumes, py::call_guard<py::gil_scoped_release>())
        .def("compute_dual_volume_k", &GeometryComplex::compute_dual_volume_k)
        .def("compute_primal_volumes", &GeometryComplex::compute_primal_volumes, py::call_guard<py::gil_scoped_release>())
        .def("compute_primal_volume_k", &GeometryComplex::compute_primal_volume_k)
        .def("compute_primal_volumes", &GeometryComplex::get_highest_dim_circumcenters, py::call_guard<py::gil_scoped_release>())
        .def("barycentric_gradients", &GeometryComplex::barycentric_gradients)
        .def("simplex_quivers", &GeometryComplex::simplex_quivers)
        .def_readonly("primal_volume", &GeometryComplex::primal_volume)
        .def_readonly("dual_volume", &GeometryComplex::dual_volume);


    py::class_<FiniteElementExteriorCalculus, GeometryComplex>(m, "FiniteElementExteriorCalculus")
		.def(py::init<>())
        .def(py::init<SimplicialComplex>())
		.def("compute_hodge_stars", &FiniteElementExteriorCalculus::compute_hodge_stars, py::call_guard<py::gil_scoped_release>())
        .def("compute_hodge_star_k", &FiniteElementExteriorCalculus::compute_hodge_star_k, py::call_guard<py::gil_scoped_release>())
		.def("set_hodge_stars_to_null", &FiniteElementExteriorCalculus::set_hodge_stars_to_null)
        .def("mass_matrix_bb_0", &FiniteElementExteriorCalculus::mass_matrix_bb_0)
        .def("bernstein", &FiniteElementExteriorCalculus::bernstein)
        .def("omega_ij", &FiniteElementExteriorCalculus::omega_ij)
        .def("compute_index_sets_o", &FiniteElementExteriorCalculus::compute_index_sets_o)
        .def("compute_index_sets_t", &FiniteElementExteriorCalculus::compute_index_sets_t)
        .def("compute_index_sets_p", &FiniteElementExteriorCalculus::compute_index_sets_p)
        .def("bb_error", &FiniteElementExteriorCalculus::bb_error, py::call_guard<py::gil_scoped_release>())
        .def_readonly("hodge_stars", &FiniteElementExteriorCalculus::hodge_stars)
        .def_readonly("all_hodge_stars", &FiniteElementExteriorCalculus::all_hodge_stars);


    py::class_<DiscreteExteriorCalculus, GeometryComplex>(m, "DiscreteExteriorCalculus")
    	.def(py::init<>())
        .def(py::init<SimplicialComplex>())
    	.def("compute_hodge_stars", &DiscreteExteriorCalculus::compute_hodge_stars)
        .def("compute_hodge_star_k", &DiscreteExteriorCalculus::compute_hodge_star_k)
    	.def("set_hodge_stars_to_null", &DiscreteExteriorCalculus::set_hodge_stars_to_null)
        .def_readonly("all_hodge_stars", &DiscreteExteriorCalculus::all_hodge_stars)
        .def_readonly("hodge_stars", &DiscreteExteriorCalculus::hodge_stars);

    m.def("get_analytical_soln", &get_analytical_soln);
    m.def("error_0", &error_0, py::call_guard<py::gil_scoped_release>());
    m.def("quadratic_error_0", &quadratic_error_0, py::call_guard<py::gil_scoped_release>());
    m.def("quadratic_error_0_bb", &quadratic_error_0_bb, py::call_guard<py::gil_scoped_release>());
    m.def("quadratic_error_0_bb_mass", &quadratic_error_0_bb, py::call_guard<py::gil_scoped_release>());
    m.def("get_simplex_volume", &get_simplex_volume);

}