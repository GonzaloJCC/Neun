/*************************************************************


*************************************************************/

#ifndef STDPSYNAPSEMODEL_H_
#define STDPSYNAPSEMODEL_H_

/**
* @brief Implements a synapse based on (Song, Miller & Abbott, 2000)
*/
template <typename precission = double>
class STDPSynapseModel {
#ifndef __AVR_ARCH__
  static_assert(std::is_floating_point<precission>::value);
#endif  //__AVR_ARCH__

  public:
    enum variable {
      g,              // Synaptic conductance (w)
      s,              // Synaptic gating variable
      n_variables
    };

    enum parameter {
      v_pre,        // Presynaptic neuron's voltage (membrane potential)
      v_post,       // Postsynaptic neuron's voltage (membrane potential)
      A_minus,      // Amplitude of synaptic modification (for LTD) (sugested 0.00525) (A- / A+ = 1.05)
      A_plus,       // Amplitude of synaptic modification (for LTP) (sugested 0.005)
      tau_minus,    // Time constant for decrease in synaptic weight (sugested 20ms)
      tau_plus,     // Time constant for increase in synaptic weight (sugested 20ms)
      g_max,        // Maximum allowed synaptic conductance
      g_min,        // Minimum allowed synaptic conductance
      i,            // Synaptic intensity
      E_syn,        // Reversal potential
      tau_syn,      // Time constant
      spike_threshold,
      n_parameters
    };
    
    typedef precission precission_t;

  public:
    STDPSynapseModel() {}

    void eval(const precission* const vars, const precission* const params,
              precission* const incs) const {

      incs[g] = 0; // d(g)/dt = 0

      incs[s] = -vars[s] / params[tau_syn]; // d(s)/dt = -s/tau_syn
      if(vars[s] < 0) incs[s] = 0;
    }
  
};

#endif /*STDPSYNAPSEMODEL_H_*/
