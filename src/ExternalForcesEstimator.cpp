#include "ExternalForcesEstimator.h"

#include <mc_control/GlobalPluginMacros.h>

namespace mc_plugin
{

ExternalForcesEstimator::~ExternalForcesEstimator() = default;

void ExternalForcesEstimator::init(mc_control::MCGlobalController & controller, const mc_rtc::Configuration & config)
{
  auto & ctl = static_cast<mc_control::MCGlobalController &>(controller);
  auto & robot = ctl.robot(ctl.robots()[0].name());
  auto & rjo = robot.refJointOrder();

  jointNumber = ctl.robot(ctl.robots()[0].name()).refJointOrder().size();

  Eigen::VectorXd qdot(jointNumber);
  for (size_t i = 0; i < jointNumber; i++)
  {
    qdot[i] = robot.alpha()[robot.jointIndexByName(rjo[i])][0];
  }

  float k = config("residual_gain", 0.0);
  referenceFrame = config("reference_frame",(std::string)"");
  Eigen::VectorXd diag = Eigen::VectorXd::Zero(jointNumber);
  diag.fill(k);

  jac = rbd::Jacobian(robot.mb(),referenceFrame);
  coriolis = new rbd::Coriolis(robot.mb());
  forwardDynamics = rbd::ForwardDynamics(robot.mb());

  residualGains = Eigen::MatrixXd(jointNumber,jointNumber);

  forwardDynamics.computeH(robot.mb(),robot.mbc());
  auto inertiaMatrix = forwardDynamics.H();
  pzero = inertiaMatrix*qdot;

  integralTerm = Eigen::VectorXd::Zero(jointNumber);
  residual = Eigen::VectorXd::Zero(jointNumber);
  externalForces = sva::ForceVecd::Zero();

  counter = 0;

  addGui(controller);
  addLog(controller);
  
  mc_rtc::log::info("ExternalForcesEstimator::init called with configuration:\n{}", config.dump(true, true));
}

void ExternalForcesEstimator::reset(mc_control::MCGlobalController & controller)
{
  removeLog(controller);
  mc_rtc::log::info("ExternalForcesEstimator::reset called");
}

void ExternalForcesEstimator::before(mc_control::MCGlobalController & controller)
{
  auto & ctl = static_cast<mc_control::MCGlobalController &>(controller);
  auto & robot = ctl.realRobot(ctl.robots()[0].name());
  auto & rjo = robot.refJointOrder();

  Eigen::VectorXd qdot(jointNumber), tau(jointNumber);
  for (size_t i = 0; i < jointNumber; i++)
  {
    qdot[i] = robot.alpha()[robot.jointIndexByName(rjo[i])][0];
    tau[i] = robot.jointTorques()[i];
  }

  forwardDynamics.computeC(robot.mb(),robot.mbc());
  forwardDynamics.computeH(robot.mb(),robot.mbc());
  auto coriolisMatrix = coriolis->coriolis(robot.mb(),robot.mbc());
  auto coriolisGravityTerm = forwardDynamics.C();
  integralTerm += (tau + (coriolisMatrix + coriolisMatrix.transpose())*qdot - coriolisGravityTerm + residual)*ctl.timestep();
  auto inertiaMatrix = forwardDynamics.H();
  auto pt = inertiaMatrix*qdot;

  residual = pt - integralTerm + pzero;
  auto pinvJtranspose = jac.jacobian(robot.mb(),robot.mbc()).transpose().completeOrthogonalDecomposition().pseudoInverse();
  externalForces = sva::ForceVecd(pinvJtranspose*residual);

  if (++counter%1000 == 10000)
  {
    mc_rtc::log::info("Tau             : {}", tau.transpose());
    mc_rtc::log::info("CoriolisGravity : {}", coriolisGravityTerm.transpose());
    mc_rtc::log::info("ExternalTorques : {}", residual.transpose());
    mc_rtc::log::info("ExternalForces  : {}", externalForces);
  }

  // mc_rtc::log::info("ExternalForcesEstimator::before");
}

void ExternalForcesEstimator::after(mc_control::MCGlobalController & controller)
{
  // mc_rtc::log::info("ExternalForcesEstimator::after");
}

mc_control::GlobalPlugin::GlobalPluginConfiguration ExternalForcesEstimator::configuration()
{
  mc_control::GlobalPlugin::GlobalPluginConfiguration out;
  out.should_run_before = true;
  out.should_run_after = false;
  out.should_always_run = true;
  return out;
}

void ExternalForcesEstimator::addGui(mc_control::MCGlobalController & controller)
{
  auto & ctl = static_cast<mc_control::MCGlobalController &>(controller);
  ctl.controller().gui()->addElement({"Plugins","External forces estimator"},
    mc_rtc::gui::Force(
      "EndEffector", [this]() { return 10*this->externalForces; },
      [this, &controller]() {
        auto transform = controller.robot().bodyPosW(referenceFrame);
        transform.rotation() = Eigen::Matrix3d::Identity();
        return transform;
      }
    )
  );
}

void ExternalForcesEstimator::addLog(mc_control::MCGlobalController & controller)
{
  controller.controller().logger().addLogEntry("ExternalForceEstimator_gain", [&, this]() {return this->residualGains(0,0);});
  controller.controller().logger().addLogEntry("ExternalForceEstimator_wrench", [&, this]() {return this->externalForces;});
  controller.controller().logger().addLogEntry("ExternalForceEstimator_residual", [&, this]() {return this->residual;});
  controller.controller().logger().addLogEntry("ExternalForceEstimator_integralTerm", [&, this]() {return this->integralTerm;});
}

void ExternalForcesEstimator::removeLog(mc_control::MCGlobalController & controller)
{
  controller.controller().logger().removeLogEntry("ExternalForceEstimator_gain");
  controller.controller().logger().removeLogEntry("ExternalForceEstimator_wrench");
  controller.controller().logger().removeLogEntry("ExternalForceEstimator_residual");
  controller.controller().logger().removeLogEntry("ExternalForceEstimator_integralTerm");
}

} // namespace mc_plugin

EXPORT_MC_RTC_PLUGIN("ExternalForcesEstimator", mc_plugin::ExternalForcesEstimator)
