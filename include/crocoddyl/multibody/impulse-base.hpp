///////////////////////////////////////////////////////////////////////////////
// BSD 3-Clause License
//
// Copyright (C) 2019-2023, LAAS-CNRS, University of Edinburgh,
//                          Heriot-Watt University
// Copyright note valid unless otherwise stated in individual files.
// All rights reserved.
///////////////////////////////////////////////////////////////////////////////

#ifndef CROCODDYL_MULTIBODY_IMPULSE_BASE_HPP_
#define CROCODDYL_MULTIBODY_IMPULSE_BASE_HPP_

#include "crocoddyl/core/utils/deprecate.hpp"
#include "crocoddyl/multibody/force-base.hpp"
#include "crocoddyl/multibody/fwd.hpp"
#include "crocoddyl/multibody/states/multibody.hpp"

namespace crocoddyl {

template <typename _Scalar>
class ImpulseModelAbstractTpl {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef _Scalar Scalar;
  typedef MathBaseTpl<Scalar> MathBase;
  typedef ImpulseDataAbstractTpl<Scalar> ImpulseDataAbstract;
  typedef StateMultibodyTpl<Scalar> StateMultibody;
  typedef typename MathBase::VectorXs VectorXs;
  typedef typename MathBase::MatrixXs MatrixXs;

  ImpulseModelAbstractTpl(boost::shared_ptr<StateMultibody> state,
                          const pinocchio::ReferenceFrame type,
                          const std::size_t nc);

  DEPRECATED(
      "Use constructor that passes the type type of contact, this assumes is "
      "pinocchio::LOCAL",
      ImpulseModelAbstractTpl(boost::shared_ptr<StateMultibody> state,
                              const std::size_t nc);)
  virtual ~ImpulseModelAbstractTpl();

  virtual void calc(const boost::shared_ptr<ImpulseDataAbstract>& data,
                    const Eigen::Ref<const VectorXs>& x) = 0;
  virtual void calcDiff(const boost::shared_ptr<ImpulseDataAbstract>& data,
                        const Eigen::Ref<const VectorXs>& x) = 0;

  virtual void updateForce(const boost::shared_ptr<ImpulseDataAbstract>& data,
                           const VectorXs& force) = 0;
  void updateForceDiff(const boost::shared_ptr<ImpulseDataAbstract>& data,
                       const MatrixXs& df_dx) const;
  void setZeroForce(const boost::shared_ptr<ImpulseDataAbstract>& data) const;
  void setZeroForceDiff(
      const boost::shared_ptr<ImpulseDataAbstract>& data) const;

  virtual boost::shared_ptr<ImpulseDataAbstract> createData(
      pinocchio::DataTpl<Scalar>* const data);

  const boost::shared_ptr<StateMultibody>& get_state() const;
  std::size_t get_nc() const;
  DEPRECATED("Use get_nc().", std::size_t get_ni() const;)
  std::size_t get_nu() const;

  std::size_t get_nf() const;

  /**
   * @brief Return the reference frame id
   */
  pinocchio::FrameIndex get_id() const;

  /**
   * @brief Modify the reference frame id
   */
  void set_id(const pinocchio::FrameIndex id);

  /**
   * @brief Modify the type of contact
   */
  void set_type(const pinocchio::ReferenceFrame type);

  /**
   * @brief Return the type of contact
   */
  pinocchio::ReferenceFrame get_type() const;

  /**
   * @brief Print information on the impulse model
   */
  template <class Scalar>
  friend std::ostream& operator<<(std::ostream& os,
                                  const ImpulseModelAbstractTpl<Scalar>& model);

  /**
   * @brief Print relevant information of the impulse model
   *
   * @param[out] os  Output stream object
   */
  virtual void print(std::ostream& os) const;

 protected:
  boost::shared_ptr<StateMultibody> state_;
  std::size_t nc_;
  pinocchio::FrameIndex id_;        //!< Reference frame id of the contact
  pinocchio::ReferenceFrame type_;  //!< Type of contact
};

template <typename _Scalar>
struct ImpulseDataAbstractTpl : public InteractionDataAbstractTpl<_Scalar> {
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef _Scalar Scalar;
  typedef MathBaseTpl<Scalar> MathBase;
  typedef typename pinocchio::DataTpl<Scalar> PinocchioData;

  typedef ForceDataAbstractTpl<Scalar> Base;
  typedef typename MathBase::VectorXs VectorXs;
  typedef typename MathBase::MatrixXs MatrixXs;
  typedef typename pinocchio::SE3Tpl<Scalar> SE3;

  template <template <typename Scalar> class Model>
  ImpulseDataAbstractTpl(Model<Scalar>* const model,
                         pinocchio::DataTpl<Scalar>* const data)
      : InteractionDataAbstractTpl<Scalar>(model, 1),
        pinocchio(data),
        dv0_dq(model->get_nc(), model->get_state()->get_nv()),
        dtau_dq(model->get_state()->get_nv(), model->get_state()->get_nv()),
        Jc(model->get_nc(), model->get_state()->get_nv()) {
    dv0_dq.setZero();
    dtau_dq.setZero();
    Jc.setZero();
  }
  virtual ~ImpulseDataAbstractTpl() {}

  PinocchioData* pinocchio;  //!< Pinocchio data

  using InteractionDataAbstractTpl<Scalar>::nf;
  using InteractionDataAbstractTpl<Scalar>::force_datas;
  using InteractionDataAbstractTpl<Scalar>::df_dx;

  MatrixXs dv0_dq;
  MatrixXs dtau_dq;

  MatrixXs Jc;  //!< Contact Jacobian
};

}  // namespace crocoddyl

/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
/* --- Details -------------------------------------------------------------- */
#include "crocoddyl/multibody/impulse-base.hxx"

#endif  // CROCODDYL_MULTIBODY_IMPULSE_BASE_HPP_
