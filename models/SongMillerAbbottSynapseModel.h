/*************************************************************


*************************************************************/

#ifndef SONGMILLERABBOTTMODEL_H_
#define SONGMILLERABBOTTMODEL_H_

/**
* @brief Implements a synapse based on (Song, Miller & Abbott, 2000)
*/
template <typename precission = double>
class SongMillerAbbottModel {
#ifndef __AVR_ARCH__
  static_assert(std::is_floating_point<precission>::value);
#endif  //__AVR_ARCH__

  public:
    enum variable {
      g,            // Synaptic conductance (w)
      n_variables
    };

    enum parameter {
      v_pre,        // Presynaptic neuron's voltage (membrane potential)
      v_post,       // Postsynaptic neuron's voltage (membrane potential)
      A_minus,      // Amplitude of synaptic modification (for LTD) (sugested 0.02)
      A_plus,       // Amplitude of synaptic modification (for LTP) (sugested 0.005)
      tau_minus,    // Time constant for decrease in synaptic weight (sugested 20ms)
      tau_plus,     // Time constant for increase in synaptic weight (sugested 20ms)
      delta_t_plus, // Time difference between presynaptic and postsynaptic spikes for increase in synaptic weight
      delta_t_minus,// Time difference between presynaptic and postsynaptic spikes for decrease in synaptic weight
      g_max,        // Maximum allowed synaptic conductance
      i,            // Synaptic intensity
      n_parameters
    };
    
    typedef precission precission_t;

  public:
    SongMillerAbbottModel() {}

    void eval(const precission* const vars, const precission* const params,
              precission* const incs) const {
        // (if delta_t > 0) {
        //   incs[g] = ...
        // } else if (delta_t < 0) {
        //   incs[g] = ...
        // } else {
        //   incs[g] = 0;
        // }
    }
  
};

#endif /*SONGMILLERABBOTTMODEL_H_*/
