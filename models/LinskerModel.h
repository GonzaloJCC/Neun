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
      v_pre,    // Presynaptic neuron's voltage (membrane potential)
      v_post,   // Postsynaptic neuron's voltage (membrane potential)
      xo,       // Constant for presynaptic neuron
      yo,       // Constant for postsynaptic neuron
      eta,      // Learning rate
      ki,       // Constant
      w_max,    // Maximum allowed synaptic weight (w_min = -wmax)
      i,        // Synaptic intensity
      n_parameters
    };

    typedef precission precission_t;

 public:
    // in constructor create mediator
    LinskerModel() {}

    void eval(const precission* const vars, const precission* const params,
              precission* const incs) const {
        // update weight

        // call mediator

        // update weight
        ;
    }
};












#endif LINSKER_MODEL_H_