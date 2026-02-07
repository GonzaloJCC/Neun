/*************************************************************


*************************************************************/

#include <DifferentialNeuronWrapper.h>
#include <SongMillerAbbottSynapse.h>
#include <HodgkinHuxleyModel.h>
#include <SystemWrapper.h>
#include <RungeKutta4.h>
#include <iostream>
#include <cstdlib>

typedef RungeKutta4 Integrator;
typedef DifferentialNeuronWrapper<SystemWrapper<HodgkinHuxleyModel<double>>, Integrator> HH;
typedef SongMillerAbbottSynapse<HH, HH, Integrator, double> Synapse;

int main(int argc, char **argv) {

  // Struct to initialize neuron model parameters
  HH::ConstructorArgs args;

  // Set the parameter values
  args.params[HH::cm] = 1 * 7.854e-3;
  args.params[HH::vna] = 50;
  args.params[HH::vk] = -77;
  args.params[HH::vl] = -54.387;
  args.params[HH::gna] = 120 * 7.854e-3;
  args.params[HH::gk] = 36 * 7.854e-3;
  args.params[HH::gl] = 0.3 * 7.854e-3;

  // Struct to initialize synapse model parameters
  Synapse::ConstructorArgs syn_args;
  syn_args.params[Synapse::A_minus] = 0.00525;
  syn_args.params[Synapse::A_plus] = 0.005;
  syn_args.params[Synapse::tau_minus] = 20;
  syn_args.params[Synapse::tau_plus] = 20;
  syn_args.params[Synapse::spike_threshold] = -54;
  syn_args.params[Synapse::g_max] = 1;
  syn_args.params[Synapse::E_syn] = 0;
  
  // syn_args.variables[Synapse::g] = 0;
  // syn_args.variables[Synapse::time_left_pre] = 0;
  // syn_args.variables[Synapse::time_left_post] = 0;


  // Set the integration step
  const double step = 0.005;
  double simulation_time = 10000;

  // Initialize neuron models
  HH h1(args), h3(args), h2(args);

  // Set initial value of V
  h1.set(HH::v, -75);
  h3.set(HH::v, -75);
  h2.set(HH::v, -75);
  
  // Initialize synapse
  Synapse s1(h1, HH::v, h2, HH::v, syn_args, 1);
  Synapse s2(h3, HH::v, h2, HH::v, syn_args, 1);

  // Initialize g
  s1.set_g(0.0052);
  s2.set_g(0.005);

  // s1.set_time_left_pre(20);
  // s1.set_time_left_post(20);
  // s2.set_time_left_pre(20);
  // s2.set_time_left_post(20);

  std::cout << "Time vpre1 vpre2 vpost i1 i2 g1 g2 SUM(g)" << std::endl;

  for (double time = 0; time < simulation_time; time += step) {
    s1.step(step, h1.get(HH::v), h2.get(HH::v));
    s2.step(step, h3.get(HH::v), h2.get(HH::v));

    h1.add_synaptic_input(0.5);
    h3.add_synaptic_input(0.5);

    h2.add_synaptic_input(0.5);
    h2.add_synaptic_input(s1.get(Synapse::i));
    h2.add_synaptic_input(s2.get(Synapse::i));
    
    h1.step(step);
    h3.step(step);
    h2.step(step);

    std::cout << time << " " << h1.get(HH::v) << " " << h3.get(HH::v) << " " << h2.get(HH::v) << " " 
              << s1.get(Synapse::i) << " " << s2.get(Synapse::i) << " " 
              << s1.get(Synapse::g) << " " << s2.get(Synapse::g) << " " 
              << s1.get(Synapse::g) + s2.get(Synapse::g) << "\n";
  }

  return 0;
}