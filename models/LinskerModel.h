/*************************************************************


*************************************************************/

#ifndef LINSKER_MODEL_H_
#define LINSKER_MODEL_H_

/**
 * @brief Implements a synapsis based on (Linsker 1986)
 */
template <typename precission = double>
class LinskerModel {
#ifndef __AVR_ARCH__
  static_assert(std::is_floating_point<precission>::value);
#endif  //__AVR_ARCH__

 public:
    // static vector<precission> all_sinapsis_i; // Possible variable for normalization
    enum variable {
      w,        // Synaptic weight (the stregth of the synapsis)
      avg_pre,  // Average of the presynaptic neuron's voltage
      avg_post, // Average of the postynaptic neuron's voltage
      avg_corr, // Average of the correlation between pre and post voltages
      n_variables
    };

    enum parameter {
      v_pre,    // Presynaptic neuron's voltage
      v_post,   // Postsynaptic neuron's voltage
      k1,       // Hebbian term coefficient
      k2,       // Presynaptic bias term coefficient
      k3,       // Postsynaptic bias term coefficient
      k4,       // Constant term coefficient (often zero)
      w_min,    // Minimum allowed synaptic weight
      w_max,    // Maximum allowed synaptic weight
      i,        // Synaptic intensity neuron
      time_constant // Constant for the averages (must have the same units as step (ms, seconds, etc.))
      n_parameters
    };

    typedef precission precission_t;

 public:
    LinskerModel() {}

    void eval(const precission* const vars, const precission* const params,
              precission* const incs) const {
        ;
    }
};












#endif LINSKER_MODEL_H_