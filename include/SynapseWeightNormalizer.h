/*************************************************************

*************************************************************/

#ifndef SYNAPSEWEIGHTNORMALIZER_H_
#define SYNAPSEWEIGHTNORMALIZER_H_

#include <algorithm>
#include <vector>
#include <cmath>
#include <map>


/**
 * @brief Implements a Singleton Mediator to normalize synapse weights
 */
template <typename Neuron, typename Synapse> // Neuron's type, synapse's type
class SynapseWeightNormalizer {

 public:

  using precission = typename Synapse::precission_t;
  static SynapseWeightNormalizer& get_instance() {
    static SynapseWeightNormalizer instance;
    return instance;
  }

  // Deny copying or reassigning the constructor (Singleton)
  SynapseWeightNormalizer(const SynapseWeightNormalizer&) = delete;
  void operator=(const SynapseWeightNormalizer&) = delete;


  /**
   * @brief Adds a synapse to the map
   */
  void add_synapse(Synapse* synapse, Neuron* post_synaptic_neuron) {
    neuron_synapses_map[post_synaptic_neuron].push_back(synapse);
  }

  /**
   * @brief Removes a synapse from the map
   */
  void remove_synapse(Synapse* synapse, Neuron* post_synaptic_neuron) {
    auto map_iterator = neuron_synapses_map.find(post_synaptic_neuron);
    if (map_iterator == neuron_synapses_map.end()) return;

    std::vector<Synapse*>& synapses = map_iterator -> second;

    auto vector_iterator = std::find(synapses.begin(), synapses.end(), synapse);

    if (vector_iterator != synapses.end()) {
      synapses.erase(vector_iterator);
    }
  }


  /**
   * @brief Upon call from a synapse which's weight has been updated normalize all synapses connected to post_synaptic_neuron
   * Synapse MUST have w_max param, get_w_max, get_weight and set_weight methods.
   */
  void normalize_weights(Synapse* updated_synapse, Neuron* post_synaptic_neuron) {
        
    auto map_iterator = neuron_synapses_map.find(post_synaptic_neuron);
    if (map_iterator == neuron_synapses_map.end()) return;

    std::vector<Synapse*>& synapse_group = map_iterator->second;
    if (synapse_group.empty()) return;

    const precission w_max = updated_synapse->get_w_max();
    const precission w_min = -w_max;
    precission range = w_max - w_min;

    // Find the real range
    precission current_max = synapse_group[0]->get_weight();
    precission current_min = synapse_group[0]->get_weight();

    for(Synapse* s : synapse_group) {
      precission w = s-> get_weight();
      if (w < current_min) current_min = w;
      if (w > current_max) current_max = w;
    }
    range = current_max - current_min;

    if (std::abs(range) < 1e-9) {
      for(Synapse* s : synapse_group) {
            precission w = s->get_weight();
            if (w > w_max) s->set_weight(w_max);
            else if (w < -w_max) s->set_weight(-w_max);
      }
      return;
    }

    // For all synapsis (except updated_synapse) do the normalization
    for (Synapse* s : synapse_group) {

        precission w = s->get_weight();
        
        // w = (w - w.min()) / (w.max() - w.min()) * 2 * wmax - wmax 
        precission new_w = (w - current_min) / (range) * 2 * w_max - w_max;
        
        s->set_weight(new_w);
    }
  }
 private:

 // Private constructor (Singleton)
  SynapseWeightNormalizer() = default;

  // Map to store the synapses connected to a neuron
  std::map<Neuron*, std::vector<Synapse*>> neuron_synapses_map;

};

#endif /*SYNAPSEWEIGHTNORMALIZER_H_*/
