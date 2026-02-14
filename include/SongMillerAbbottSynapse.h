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

#define CURRENT_DIRECTION -1

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
        System(synapse) {
          for(int i = 0; i < System::n_variables; i++) {
            System::m_variables[i] = 0;
          }
        }
 private:
  void calculate_i() {
    precission E_syn = System::m_parameters[System::E_syn];
  
  // Isyn = gsyn * s * (V - Esyn)

  System::m_parameters[System::i] = CURRENT_DIRECTION * System::m_variables[System::g] * System::m_variables[System::s] * (System::m_parameters[System::v_post] - E_syn);
  }

  void update_g(precission h) {

    precission threshold = System::m_parameters[System::spike_threshold];
    int activation_flag = 0; // Flag to update g only if the spike has just occurred
    // Detect spikes // ----------------------------------------------------------------------------------
     
    // For vpre
    if (System::m_parameters[System::v_pre] >= threshold && System::m_variables[System::time_left_pre] <= 0) {
      System::m_variables[System::time_left_pre] = System::m_parameters[System::tau_plus];
      activation_flag = 1;
      System::m_variables[System::s] = 1.0; // Update s if vpre spikes
    }

    // For vpost
    if (System::m_parameters[System::v_post] >= threshold && System::m_variables[System::time_left_post] <= 0) {
      System::m_variables[System::time_left_post] = System::m_parameters[System::tau_minus];
      activation_flag = 1;
    }
    // ---------------------------------------------------------------------------------------------------

    precission delta_t = 0;
    precission f_delta_t = 0;
    // When both neurons spiked
    if (activation_flag && System::m_variables[System::time_left_post] > 0 && System::m_variables[System::time_left_pre] > 0) {
      delta_t = System::m_variables[System::time_left_pre] - System::m_variables[System::time_left_post];
      // △t = t_pre - t_post
    }

    // LTP {F(△t) = (A+ * exp(△t/ τ+)) if △t < 0}
    if (delta_t < 0) {
      f_delta_t = (System::m_parameters[System::A_plus] * std::exp(delta_t / System::m_parameters[System::tau_plus]));
    }

    // LTD {F(△t) = (-A- * exp(-△t/ τ-)) if △t > 0}
    if (delta_t > 0) {
      f_delta_t = (-System::m_parameters[System::A_minus] * std::exp(-delta_t / System::m_parameters[System::tau_minus]));
    }

    System::m_variables[System::g] += f_delta_t; // g := g + F(△t)

    /* Limit growth */
    if (System::m_variables[System::g] > System::m_parameters[System::g_max]) {
      System::m_variables[System::g] = System::m_parameters[System::g_max];
    }

    /* Limit decay */
    if (System::m_variables[System::g] < System::m_parameters[System::g_min]) {
      System::m_variables[System::g] = System::m_parameters[System::g_min];
    }
  }

 public:
  void step(precission h) {
    
    precission v_pre = m_n1.get(m_n1_variable);
    System::m_parameters[System::v_pre] = v_pre;
    precission v_post = m_n2.get(m_n2_variable);
    System::m_parameters[System::v_post] = v_post;

    for (int i = 0; i < m_steps; ++i) {
      TIntegrator::step(*this, h, System::m_variables, System::m_parameters);
    }

    update_g(h);

    // Calculate synaptic current
    calculate_i();

  }

  void step(precission h, precission vpre, precission vpost) {
    
    System::m_parameters[System::v_pre] = vpre;
    System::m_parameters[System::v_post] = vpost;

    for (int i = 0; i < m_steps; ++i) {
      TIntegrator::step(*this, h, System::m_variables, System::m_parameters);
    }

    update_g(h);

    // Calculate synaptic current
    calculate_i();
  }

  void set_g(precission g) {
    System::m_variables[System::g] = g;
  }

  void set_time_left_pre(precission time_left_pre) {
    System::m_variables[System::time_left_pre] = time_left_pre;
  }

  void set_time_left_post(precission time_left_post) {
    System::m_variables[System::time_left_post] = time_left_post;
  }
};

#endif /*SONGMILLERABBOTTSYNAPSE_H_*/
