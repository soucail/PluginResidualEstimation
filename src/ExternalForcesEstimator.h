/*
 * Copyright 2021 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

#include <mc_control/GlobalPlugin.h>
#include <mc_rtc/log/FlatLog.h>

#include <RBDyn/Coriolis.h>
#include <RBDyn/MultiBody.h>
#include <RBDyn/MultiBodyConfig.h>

#include <mc_rbdyn/ExternalTorqueSensor.h>
#include <mc_rbdyn/VirtualTorqueSensor.h>
#include <mc_rbdyn/Robot.h>


#include "utils/ROSSubscriber.h"

namespace mc_plugin
{

struct ExternalForcesEstimator : public mc_control::GlobalPlugin
{
  void init(mc_control::MCGlobalController & controller, const mc_rtc::Configuration & config) override;

  void reset(mc_control::MCGlobalController & controller) override;

  void before(mc_control::MCGlobalController &) override;

  void after(mc_control::MCGlobalController & controller) override;

  mc_control::GlobalPlugin::GlobalPluginConfiguration configuration() override;

  ~ExternalForcesEstimator() override;

  void rosSpinner(void);

  void addGui(mc_control::MCGlobalController & controller);
  void addLog(mc_control::MCGlobalController & controller);
  void removeLog(mc_control::MCGlobalController & controller);

private:
  int jointNumber;
  int counter;
  double dt;
  bool verbose;
  bool isActive;

  double residualGains;
  double power_;

  std::string referenceFrame;

  rbd::Jacobian jac;
  rbd::Coriolis * coriolis;
  rbd::ForwardDynamics forwardDynamics;

  Eigen::VectorXd pzero;

  Eigen::VectorXd integralTerm;
  Eigen::VectorXd residual;
  Eigen::VectorXd FTSensorTorques;
  Eigen::VectorXd prevFTSensorTorques;
  Eigen::VectorXd filteredFTSensorTorques;
  Eigen::VectorXd newExternalTorques;
  Eigen::VectorXd externalTorques;
  Eigen::VectorXd filteredExternalTorques;
  sva::ForceVecd externalForces;
  sva::ForceVecd externalForcesResidual;
  sva::ForceVecd newExternalForces;
  sva::ForceVecd filteredFTSensorForces;
  Eigen::Vector6d externalForcesFT;
  mc_rbdyn::ExternalTorqueSensor * extTorqueSensor;
  mc_rbdyn::VirtualTorqueSensor * virtTorqueSensor;
  Eigen::VectorXd motor_current_;
  Eigen::VectorXd tau_current_;
  Eigen::VectorXd alpha_;


  // Force sensor
  bool use_force_sensor_;

  bool ros_force_sensor_;
  std::shared_ptr<ros::NodeHandle> nh_;
  std::thread spinThread_;
  double maxTime_ = 0.001;
  double freq_ = 1000;
  std::string force_sensor_topic_ = "/fast_chatter";
  ROSWrenchStampedSubscriber wrench_sub_;
};

} // namespace mc_plugin
