/*************************************************************


*************************************************************/

#include <DifferentialNeuronWrapper.h>
#include <LinskerSynapsis.h>
#include <HodgkinHuxleyModel.h>
#include <SystemWrapper.h>
#include <RungeKutta4.h>
#include <iostream>
#include <cstdlib>

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


  // Synapsis::ConstructorArgs syn_args;
  // syn_args.params[Synapsis::xo] = -65;
  // syn_args.params[Synapsis::yo] = -65;
  // syn_args.params[Synapsis::eta] = 0.00001;
  // syn_args.params[Synapsis::k1] = -500;
  // syn_args.params[Synapsis::w_max] = 3;

  Synapsis::ConstructorArgs syn_args;
  syn_args.params[Synapsis::xo] = -65;
  syn_args.params[Synapsis::yo] = -63;
  syn_args.params[Synapsis::eta] = 0.00001;
  syn_args.params[Synapsis::k1] = -50;
  syn_args.params[Synapsis::w_max] = 2;


  // Set the integration step
  const double step = 0.005;
  double simulation_time = 10000;

  // 3 Neurons 2 synapses
  {
      // Initialize neuron models
      HH h1(args), h2(args);
      HH h3(args);

      // Set initial value of V
      h1.set(HH::v, -75);
      h3.set(HH::v, -71);
      
      // Initialize 2 synapsis
      Synapsis s1(h1, HH::v, h2, HH::v, syn_args, 1);
      Synapsis s2(h3, HH::v, h2, HH::v, syn_args, 1);

      // Initialize weights
      s1.set_weight(0.001 * (rand() / (double)RAND_MAX));
      s2.set_weight(0.001 * (rand() / (double)RAND_MAX));

      std::cout << "Time v1pre v2pre vpost i1 i2 w1 w2" << std::endl;

      for (double time = 0; time < simulation_time; time += step) {
          s1.step(step, h1.get(HH::v), h2.get(HH::v));
          s2.step(step, h3.get(HH::v), h2.get(HH::v));

          // Inputs
          h1.add_synaptic_input(0.5);
          h2.add_synaptic_input(0.5);
          h3.add_synaptic_input(0.6);

          h2.add_synaptic_input(s1.get(Synapsis::i));
          h2.add_synaptic_input(s2.get(Synapsis::i));

          h1.step(step);
          h3.step(step);
          h2.step(step);

          std::cout << time << " " << h1.get(HH::v) << " " << h3.get(HH::v) << " " << h2.get(HH::v) << " " 
                    << s1.get(Synapsis::i) << " " << s2.get(Synapsis::i) << " " 
                    << s1.get(Synapsis::w) << " " << s2.get(Synapsis::w) << " " 
                    << (s1.get(Synapsis::w) + s2.get(Synapsis::w)) << "\n";
      }
  }

  // ----------------------------------------------- //
  
  // 5 Neurons 4 synapses
//   {
//       // Initialize neuron models
//       HH h1(args), h2(args), h3(args), h4(args), h5(args);

//       // Set initial value of v
//       h1.set(HH::v, -75);
//       h3.set(HH::v, -71);
//       h4.set(HH::v, -70);
//       h5.set(HH::v, -78);

//       // Initialize 4 synapses (h2 is post-synaptic)
//       Synapsis s1(h1, HH::v, h2, HH::v, syn_args, 1);
//       Synapsis s2(h3, HH::v, h2, HH::v, syn_args, 1);
//       Synapsis s3(h4, HH::v, h2, HH::v, syn_args, 1);
//       Synapsis s4(h5, HH::v, h2, HH::v, syn_args, 1);

//       // Initialize weights
//       s1.set_weight(0.001 * (rand() / (double)RAND_MAX));
//       s2.set_weight(0.001 * (rand() / (double)RAND_MAX));
//       s3.set_weight(0.001 * (rand() / (double)RAND_MAX));
//       s4.set_weight(0.001 * (rand() / (double)RAND_MAX));

//       std::cout << "Time V1pre V2pre V3pre V4pre Vpost i1 i2 i3 i4 w1 w2 w3 w4 SUM(W)" << std::endl;

//       for (double time = 0; time < simulation_time; time += step) {
//           s1.step(step, h1.get(HH::v), h2.get(HH::v));
//           s2.step(step, h3.get(HH::v), h2.get(HH::v));
//           s3.step(step, h4.get(HH::v), h2.get(HH::v));
//           s4.step(step, h5.get(HH::v), h2.get(HH::v));

//           // External Inputs
//           h1.add_synaptic_input(0.5);
//           h2.add_synaptic_input(0.5);
//           h3.add_synaptic_input(0.6);
//           h4.add_synaptic_input(0.5);
//           h5.add_synaptic_input(0.5);

//           // Synaptic inputs to h2
//           h2.add_synaptic_input(s1.get(Synapsis::i));
//           h2.add_synaptic_input(s2.get(Synapsis::i));
//           h2.add_synaptic_input(s3.get(Synapsis::i));
//           h2.add_synaptic_input(s4.get(Synapsis::i));

//           h1.step(step);
//           h3.step(step);
//           h4.step(step);
//           h5.step(step);
//           h2.step(step);

//           std::cout << time << " " 
//                     << h1.get(HH::v) << " " << h3.get(HH::v) << " " << h4.get(HH::v) << " " << h5.get(HH::v) << " " << h2.get(HH::v) << " " 
//                     << s1.get(Synapsis::i) << " " << s2.get(Synapsis::i) << " " << s3.get(Synapsis::i) << " " << s4.get(Synapsis::i) << " " 
//                     << s1.get(Synapsis::w) << " " << s2.get(Synapsis::w) << " " << s3.get(Synapsis::w) << " " << s4.get(Synapsis::w) << " "
//                     << (s1.get(Synapsis::w) + s2.get(Synapsis::w) + s3.get(Synapsis::w) + s4.get(Synapsis::w))
//                     << "\n";
//       }
//   }

  return 0;
}
