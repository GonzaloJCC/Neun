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
      w,        // Synaptic efficacy
      n_variables
    };

    enum parameter {
      v_pre,    // Presynaptic neuron's voltage
      v_post,   // Postsynaptic neuron's voltage
      xo,       // Constant for presynaptic neuron
      yo,       // Constant for postsynaptic neuron
      eta,      // Learning rate
      ki,       // Constant
      w_min,    // Minimum allowed synaptic weight
      w_max,    // Maximum allowed synaptic weight
      i,        // Synaptic intensity neuron
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