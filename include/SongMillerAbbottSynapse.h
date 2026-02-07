/*************************************************************

*************************************************************/

#ifndef SONGMILLERABBOTTSYNAPSE_H_
#define SONGMILLERABBOTTSYNAPSE_H_

#ifndef __AVR_ARCH__
#include <type_traits>
#include <concepts>
#include "IntegratorConcept.h"
#include "NeuronConcept.h"
#endif //__AVR_ARCH__

#include "SongMillerAbbottSynapseModel.h"
#include "IntegratedSystemWrapper.h"
#include "SerializableWrapper.h"
#include "SystemWrapper.h"
#include <cmath>


/**
* Implements a synapse based on (Song, Miller & Abbott, 2000)
*/

template <typename TNode1, typename TNode2, typename TIntegrator, typename precission = double>
requires NeuronConcept<TNode1> && NeuronConcept<TNode2> &&
    IntegratorConcept<TIntegrator, SerializableWrapper<
          SystemWrapper<SongMillerAbbottSynapseModel<precission> > > >

class SongMillerAbbottSynapse : public SerializableWrapper<
          SystemWrapper<SongMillerAbbottSynapseModel<precission> > > {
 private:
  #ifndef __AVR_ARCH__
    static_assert(std::is_floating_point<precission>::value);
  #endif  //__AVR_ARCH__

  precission m_release_time;

  TNode1 const &m_n1;
  TNode2 &m_n2;

  precission m_last_value_pre;
  precission m_last_value_post;

  const typename TNode1::variable m_n1_variable;
  const typename TNode2::variable m_n2_variable;


  typedef SerializableWrapper<
      SystemWrapper<SongMillerAbbottSynapseModel<precission> > > System;

  const int m_steps;

 public:
  typedef typename System::precission_t precission_t;
  typedef typename System::variable variable;
  typedef typename System::parameter parameter;
  typedef typename System::ConstructorArgs ConstructorArgs;

  SongMillerAbbottSynapse (TNode1 const &n1, typename TNode1::variable v1,
                                              TNode2 &n2, typename TNode2::variable v2, 
                                              ConstructorArgs &args, int steps)
      : m_n1(n1),
        m_n2(n2),
        m_n1_variable(v1),
        m_n2_variable(v2),
        System(args),
        m_last_value_pre(-75),
        m_last_value_post(-75),
        m_steps(steps) {
          for(int i = 0; i < System::n_variables; i++) {
            System::m_variables[i] = 0;
          }
        }

  SongMillerAbbottSynapse (TNode1 const &n1, TNode2 &n2,
                            SongMillerAbbottSynapse  const &synapse)
      : m_n1(n1),
        m_n2(n2),
        m_n1_variable(synapse.m_n1_variable),
        m_n2_variable(synapse.m_n2_variable),
        m_steps(synapse.m_steps),
        m_last_value_pre(synapse.m_last_value_pre),
        m_last_value_post(synapse.m_last_value_post),
        System(synapse) {
          for(int i = 0; i < System::n_variables; i++) {
            System::m_variables[i] = 0;
          }
        }

  void step(precission h) {
    
    precission v_pre = m_n1.get(m_n1_variable);
    System::m_parameters[System::v_pre] = v_pre;
    precission v_post = m_n2.get(m_n2_variable);
    System::m_parameters[System::v_post] = v_post;
    
    System::m_parameters[System::i] = System::m_variables[System::g] * System::m_parameters[System::v_pre];
    precission threshold = System::m_parameters[System::spike_threshold];

    for (int i = 0; i < m_steps; ++i) {
      TIntegrator::step(*this, h, System::m_variables, System::m_parameters);
    }

    // LTP (if vpost spikes and vpre was active)
    if (m_last_value_post < threshold && v_post >= threshold) {
      precission time_left = System::m_variables[System::time_left_pre];
      if (time_left > 0) {
        System::m_variables[System::g] += System::m_parameters[System::g_max] * (System::m_parameters[System::A_plus] * std::exp(time_left / System::m_parameters[System::tau_plus]));
      } // g := g + g_max * (A+ * exp(△t/ τ+))

      // Restart time left until forgetting post-synaptic spike
      System::m_variables[System::time_left_post] = System::m_parameters[System::tau_minus];
    }

    // LTD (if vpre spikes and vpost was active)
    if (m_last_value_pre < threshold && v_pre >= threshold) {
      precission time_left = System::m_variables[System::time_left_post];
      if (time_left > 0) {
        System::m_variables[System::g] -= System::m_parameters[System::g_max] * (System::m_parameters[System::A_minus] * std::exp(time_left / System::m_parameters[System::tau_minus]));
      } // g := g - g_max * (A- * exp(△t/ τ-))

      // Restart time left until forgetting pre-synaptic spike
      System::m_variables[System::time_left_pre] = System::m_parameters[System::tau_plus];
    }

    /* Limit growth */
    if (System::m_variables[System::g] > System::m_parameters[System::g_max]) {
      System::m_variables[System::g] = System::m_parameters[System::g_max];
    }

    /* Limit decay */
    if (System::m_variables[System::g] < 0) {
      System::m_variables[System::g] = 0;
    }

    // Save states
    m_last_value_pre = v_pre;
    m_last_value_post = v_post;

  }

  void step(precission h, precission vpre, precission vpost) {
    
    System::m_parameters[System::v_pre] = vpre;
    System::m_parameters[System::v_post] = vpost;

    System::m_parameters[System::i] = System::m_variables[System::g] * System::m_parameters[System::v_pre];

    precission threshold = System::m_parameters[System::spike_threshold];

    for (int i = 0; i < m_steps; ++i) {
      TIntegrator::step(*this, h, System::m_variables, System::m_parameters);
    }

    // LTP (if vpost spikes and vpre was active)
    if (m_last_value_post < threshold && vpost >= threshold) {
      precission time_left = System::m_variables[System::time_left_pre];
      if (time_left > 0) {
        System::m_variables[System::g] += System::m_parameters[System::g_max] * (System::m_parameters[System::A_plus] * std::exp(time_left / System::m_parameters[System::tau_plus]));
      } // g := g + g_max * (A+ * exp(△t/ τ+))

      // Restart time left until forgetting post-synaptic spike
      System::m_variables[System::time_left_post] = System::m_parameters[System::tau_minus];
    }

    // LTD (if vpre spikes and vpost was active)
    if (m_last_value_pre < threshold && vpre >= threshold) {
      precission time_left = System::m_variables[System::time_left_post];
      if (time_left > 0) {
        System::m_variables[System::g] -= System::m_parameters[System::g_max] * (System::m_parameters[System::A_minus] * std::exp(time_left / System::m_parameters[System::tau_minus]));
      } // g := g - g_max * (A- * exp(△t/ τ-))

      // Restart time left until forgetting pre-synaptic spike
      System::m_variables[System::time_left_pre] = System::m_parameters[System::tau_plus];
    }

    /* Limit growth */
    if (System::m_variables[System::g] > System::m_parameters[System::g_max]) {
      System::m_variables[System::g] = System::m_parameters[System::g_max];
    }

    /* Limit decay */
    if (System::m_variables[System::g] < 0) {
      System::m_variables[System::g] = 0;
    }

    // Save states for the next step
    m_last_value_pre = vpre;
    m_last_value_post = vpost;

  }
};

#endif /*SONGMILLERABBOTTSYNAPSE_H_*/
