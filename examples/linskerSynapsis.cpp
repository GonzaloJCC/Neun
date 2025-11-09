/*************************************************************


*************************************************************/

#include <DifferentialNeuronWrapper.h>
#include <LinskerSynapsis.h>
#include <HodgkinHuxleyModel.h>
#include <SystemWrapper.h>
#include <RungeKutta4.h>
#include <iostream>

typedef RungeKutta4 Integrator;
typedef DifferentialNeuronWrapper<SystemWrapper<HodgkinHuxleyModel<double>>, Integrator> HH;
typedef LinskerSynapsis<HH, HH, Integrator, double> Synapsis;
// typedef LinskerSynapsisModel<double> SynapsisModel;

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


  Synapsis::ConstructorArgs syn_args;
  syn_args.params[Synapsis::xo] = 0.2;
  syn_args.params[Synapsis::yo] = -50;
  syn_args.params[Synapsis::eta] = 0.0005;
  syn_args.params[Synapsis::k1] = 1;
  syn_args.params[Synapsis::w_max] = 3;


  // Initialize neuron models
  HH h1(args), h2(args);

  // Set initial value of V in neuron n1
  h1.set(HH::v, -75);

  // Set the integration step
  const double step = 0.005;

  // Initialize a synapsis between the neurons
  Synapsis s(h1, HH::v, h2, HH::v, syn_args, 1);


  // Perform the simulation
  double simulation_time = 10000;
  std::cout << "Time" << " " << "Vpre" << " " << "Vpost" 
              << " " << "i" << " " << "w" << std::endl;

  for (double time = 0; time < simulation_time; time += step) {
    s.step(step, h1.get(HH::v), h2.get(HH::v));

    // Provide an external current input to both neurons
    h1.add_synaptic_input(0.5);
    h2.add_synaptic_input(0.5);

    h2.add_synaptic_input(s.get(Synapsis::i));

    h1.step(step);
    h2.step(step);

    std::cout << time << " " << h1.get(HH::v) << " " << h2.get(HH::v) 
              << " " << s.get(Synapsis::i)<< " " << s.get(Synapsis::w) << "\n";
  }

  return 0;
}

/* TODO: Add h3 and synapse s2 from h3 to h2*/