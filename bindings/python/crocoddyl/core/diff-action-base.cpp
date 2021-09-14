///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2019-2021, LAAS-CNRS, University of Edinburgh, University of Oxford
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#include "python/crocoddyl/core/core.hpp"
#include "python/crocoddyl/core/diff-action-base.hpp"
#include "python/crocoddyl/utils/printable.hpp"
#include "python/crocoddyl/utils/vector-converter.hpp"

namespace crocoddyl {
namespace python {

void exposeDifferentialActionAbstract() {
  // Register custom converters between std::vector and Python list
  typedef boost::shared_ptr<DifferentialActionModelAbstract> DifferentialActionModelPtr;
  typedef boost::shared_ptr<DifferentialActionDataAbstract> DifferentialActionDataPtr;
  StdVectorPythonVisitor<DifferentialActionModelPtr, std::allocator<DifferentialActionModelPtr>, true>::expose(
      "StdVec_DiffActionModel");
  StdVectorPythonVisitor<DifferentialActionDataPtr, std::allocator<DifferentialActionDataPtr>, true>::expose(
      "StdVec_DiffActionData");

  bp::register_ptr_to_python<boost::shared_ptr<DifferentialActionModelAbstract> >();

  bp::class_<DifferentialActionModelAbstract_wrap, boost::noncopyable>(
      "DifferentialActionModelAbstract",
      "Abstract class for the differential action model.\n\n"
      "A differential action model is the time-continuous version of an action model. Each\n"
      "node, in our optimal control problem, is described through an action model. Every\n"
      "time that we want describe a problem, we need to provide ways of computing the\n"
      "dynamics, cost functions and their derivatives. These computations are mainly carrying\n"
      "on inside calc() and calcDiff(), respectively.",
      bp::init<boost::shared_ptr<StateAbstract>, int, bp::optional<int> >(
          bp::args("self", "state", "nu", "nr"),
          "Initialize the differential action model.\n\n"
          "You can also describe autonomous systems by setting nu = 0.\n"
          ":param state: state\n"
          ":param nu: dimension of control vector\n"
          ":param nr: dimension of cost-residual vector (default 1)"))
      .def("calc", pure_virtual(&DifferentialActionModelAbstract_wrap::calc), bp::args("self", "data", "x", "u"),
           "Compute the system acceleration and cost value.\n\n"
           ":param data: differential action data\n"
           ":param x: state point (dim. state.nx)\n"
           ":param u: control input (dim. nu)")
      .def<void (DifferentialActionModelAbstract::*)(const boost::shared_ptr<DifferentialActionDataAbstract>&,
                                                     const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calc", &DifferentialActionModelAbstract::calc, bp::args("self", "data", "x"),
          "Compute the total cost value for nodes that depends only on the state.\n\n"
          "It updates the total cost and the system acceleration is not updated as it is expected to be zero.\n"
          "This function is commonly used in the terminal nodes of an optimal control problem.\n"
          ":param data: differential action data\n"
          ":param x: state point (dim. state.nx)")
      .def("calcDiff", pure_virtual(&DifferentialActionModelAbstract_wrap::calcDiff),
           bp::args("self", "data", "x", "u"),
           "Compute the derivatives of the dynamics and cost functions.\n\n"
           "It computes the partial derivatives of the dynamical system and the cost\n"
           "function. It assumes that calc has been run first.\n"
           "This function builds a quadratic approximation of the\n"
           "time-continuous action model (i.e. dynamical system and cost function).\n"
           ":param data: differential action data\n"
           ":param x: state point (dim. state.nx)\n"
           ":param u: control input (dim. nu)")
      .def<void (DifferentialActionModelAbstract::*)(const boost::shared_ptr<DifferentialActionDataAbstract>&,
                                                     const Eigen::Ref<const Eigen::VectorXd>&)>(
          "calcDiff", &DifferentialActionModelAbstract::calcDiff, bp::args("self", "data", "x"),
          "Compute the derivatives of the cost functions with respect to the state only.\n\n"
          "It updates the derivatives of the cost function with respect to the state only.\n"
          "This function is commonly used in the terminal nodes of an optimal control problem.\n"
          ":param data: action data\n"
          ":param x: state point (dim. state.nx)")
      .def("createData", &DifferentialActionModelAbstract_wrap::createData,
           &DifferentialActionModelAbstract_wrap::default_createData, bp::args("self"),
           "Create the differential action data.\n\n"
           "Each differential action model has its own data that needs to be\n"
           "allocated. This function returns the allocated data for a predefined\n"
           "DAM.\n"
           ":return DAM data.")
      .def("quasiStatic", &DifferentialActionModelAbstract_wrap::quasiStatic_x,
           DifferentialActionModel_quasiStatic_wraps(
               bp::args("self", "data", "x", "maxiter", "tol"),
               "Compute the quasic-static control given a state.\n\n"
               "It runs an iterative Newton step in order to compute the quasic-static regime\n"
               "given a state configuration.\n"
               ":param data: action data\n"
               ":param x: discrete-time state vector\n"
               ":param maxiter: maximum allowed number of iterations\n"
               ":param tol: stopping tolerance criteria (default 1e-9)\n"
               ":return u: quasic-static control"))
      .def("quasiStatic", &DifferentialActionModelAbstract_wrap::quasiStatic,
           &DifferentialActionModelAbstract_wrap::default_quasiStatic,
           bp::args("self", "data", "u", "x", "maxiter", "tol"))
      .add_property("nu", bp::make_function(&DifferentialActionModelAbstract_wrap::get_nu),
                    "dimension of control vector")
      .add_property("nr", bp::make_function(&DifferentialActionModelAbstract_wrap::get_nr),
                    "dimension of cost-residual vector")
      .add_property("state",
                    bp::make_function(&DifferentialActionModelAbstract_wrap::get_state,
                                      bp::return_value_policy<bp::return_by_value>()),
                    "state")
      .add_property("has_control_limits",
                    bp::make_function(&DifferentialActionModelAbstract_wrap::get_has_control_limits),
                    "indicates whether problem has finite control limits")
      .add_property("u_lb",
                    bp::make_function(&DifferentialActionModelAbstract_wrap::get_u_lb,
                                      bp::return_value_policy<bp::return_by_value>()),
                    &DifferentialActionModelAbstract_wrap::set_u_lb, "lower control limits")
      .add_property("u_ub",
                    bp::make_function(&DifferentialActionModelAbstract_wrap::get_u_ub,
                                      bp::return_value_policy<bp::return_by_value>()),
                    &DifferentialActionModelAbstract_wrap::set_u_ub, "upper control limits")
      .def(PrintableVisitor<DifferentialActionModelAbstract>());

  bp::register_ptr_to_python<boost::shared_ptr<DifferentialActionDataAbstract> >();

  bp::class_<DifferentialActionDataAbstract, boost::noncopyable>(
      "DifferentialActionDataAbstract",
      "Abstract class for differential action data.\n\n"
      "In crocoddyl, an action data contains all the required information for processing an\n"
      "user-defined action model. The action data typically is allocated onces by running\n"
      "model.createData() and contains the first- and second- order derivatives of the dynamics\n"
      "and cost function, respectively.",
      bp::init<DifferentialActionModelAbstract*>(
          bp::args("self", "model"),
          "Create common data shared between DAMs.\n\n"
          "The differential action data uses the model in order to first process it.\n"
          ":param model: differential action model"))
      .add_property(
          "cost",
          bp::make_getter(&DifferentialActionDataAbstract::cost, bp::return_value_policy<bp::return_by_value>()),
          bp::make_setter(&DifferentialActionDataAbstract::cost), "cost value")
      .add_property("xout", bp::make_getter(&DifferentialActionDataAbstract::xout, bp::return_internal_reference<>()),
                    bp::make_setter(&DifferentialActionDataAbstract::xout), "evolution state")
      .add_property("r", bp::make_getter(&DifferentialActionDataAbstract::r, bp::return_internal_reference<>()),
                    bp::make_setter(&DifferentialActionDataAbstract::r), "cost residual")
      .add_property("Fx", bp::make_getter(&DifferentialActionDataAbstract::Fx, bp::return_internal_reference<>()),
                    bp::make_setter(&DifferentialActionDataAbstract::Fx), "Jacobian of the dynamics")
      .add_property("Fu", bp::make_getter(&DifferentialActionDataAbstract::Fu, bp::return_internal_reference<>()),
                    bp::make_setter(&DifferentialActionDataAbstract::Fu), "Jacobian of the dynamics")
      .add_property("Lx", bp::make_getter(&DifferentialActionDataAbstract::Lx, bp::return_internal_reference<>()),
                    bp::make_setter(&DifferentialActionDataAbstract::Lx), "Jacobian of the cost")
      .add_property("Lu", bp::make_getter(&DifferentialActionDataAbstract::Lu, bp::return_internal_reference<>()),
                    bp::make_setter(&DifferentialActionDataAbstract::Lu), "Jacobian of the cost")
      .add_property("Lxx", bp::make_getter(&DifferentialActionDataAbstract::Lxx, bp::return_internal_reference<>()),
                    bp::make_setter(&DifferentialActionDataAbstract::Lxx), "Hessian of the cost")
      .add_property("Lxu", bp::make_getter(&DifferentialActionDataAbstract::Lxu, bp::return_internal_reference<>()),
                    bp::make_setter(&DifferentialActionDataAbstract::Lxu), "Hessian of the cost")
      .add_property("Luu", bp::make_getter(&DifferentialActionDataAbstract::Luu, bp::return_internal_reference<>()),
                    bp::make_setter(&DifferentialActionDataAbstract::Luu), "Hessian of the cost");
}

}  // namespace python
}  // namespace crocoddyl
