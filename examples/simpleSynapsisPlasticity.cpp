/*************************************************************


*************************************************************/

#include <DifferentialNeuronWrapper.h>
#include <SimpleSynapsisPlasticity.h>
#include <HodgkinHuxleyModel.h>
#include <SystemWrapper.h>
#include <RungeKutta4.h> 

typedef RungeKutta4 Integrator;
typedef HodgkinHuxleyModel<double> HHModel;
typedef DifferentialNeuronWrapper<SystemWrapper<HHModel>, Integrator> HH;
typedef SimpleSynapsisPlasticity<HH, HH, Integrator, double> Synapsis;

int main(int argc, char **argv) {

  // Struct to initialize neuron model parameters
  HH::ConstructorArgs args;

  // Set the parameter values for the neurons
  args.params[HH::cm] = 1 * 7.854e-3;
  args.params[HH::vna] = 50;
  args.params[HH::vk] = -77;
  args.params[HH::vl] = -54.387;
  args.params[HH::gna] = 120 * 7.854e-3;
  args.params[HH::gk] = 36 * 7.854e-3;
  args.params[HH::gl] = 0.3 * 7.854e-3;

  // Set the parameter values for the synapsis
  Synapsis::ConstructorArgs syn_args;
  syn_args.params[Synapsis::gfast] = 0.015;
  syn_args.params[Synapsis::Esyn] = -75;
  syn_args.params[Synapsis::sfast] = 0.2;
  syn_args.params[Synapsis::Vfast] = -50;
  syn_args.params[Synapsis::steps_for_degradation] = 5;
  syn_args.params[Synapsis::steps_remaining] = 5;
  syn_args.params[Synapsis::previous_vpre] = 0;

  // Initialize neuron models
  HH h1(args), h2(args);

  // Set initial value of V in neruron 1
  h1.set(HH::v, -75);

  // Set the step
  const double step = 0.005;

  // Initialize the synapsis
  // Neuron1, Neuron1 v, Neuron2, Neuron2 v, Synapsis args, steps
  Synapsis s(h1, HH::v, h2, HH::v, syn_args, 1);

  // Perform the simulation
  double simulation_time = 2000;
  std::cout << "Time" << " " << "Vpre" << " " << "Vpost" 
              << " " << "i" << " " << "gfast"
              << std::endl;

  for (double time = 0; time < simulation_time; time += step) {
    s.step(step, h1.get(HH::v), h2.get(HH::v));

    // Provide an external current input to both neurons
    h1.add_synaptic_input(0.5);
    h2.add_synaptic_input(0.5);

    h2.add_synaptic_input(s.get(Synapsis::i));

    h1.step(step);
    h2.step(step);

    std::cout << time << " " << h1.get(HH::v) << " " << h2.get(HH::v)
              << " " << s.get(Synapsis::i) << " " << s.get(Synapsis::gfast) << "\n";

  }

  return 0;
}
