/*************************************************************


*************************************************************/

#ifndef SONGMILLERABBOTTSYNAPSEMODEL_H_
#define SONGMILLERABBOTTSYNAPSEMODEL_H_

/**
* @brief Implements a synapse based on (Song, Miller & Abbott, 2000)
*/
template <typename precission = double>
class SongMillerAbbottSynapseModel {
#ifndef __AVR_ARCH__
  static_assert(std::is_floating_point<precission>::value);
#endif  //__AVR_ARCH__

  public:
    enum variable {
      g,              // Synaptic conductance (w)
      time_left_pre,  // Time left until synapse forgets presynaptic spike
      time_left_post, // Time left until synapse forgets postsynaptic spike
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

      spike_threshold,
      n_parameters
    };
    
    typedef precission precission_t;

  public:
    SongMillerAbbottSynapseModel() {}

    /* Decrement the time left until the synapse forgets the presynaptic or postsynaptic spike*/
    void eval(const precission* const vars, const precission* const params,
              precission* const incs) const {

      if (vars[time_left_pre] > 0) incs[time_left_pre] = -1; // d(time_left_pre)/dt = -1
      else incs[time_left_pre] = 0;

      if (vars[time_left_post] > 0) incs[time_left_post] = -1; // d(time_left_post)/dt = -1
      else incs[time_left_post] = 0;

      incs[g] = 0; // d(g)/dt = 0
    }
  
};

#endif /*SONGMILLERABBOTTSYNAPSEMODEL_H_*/
