/*************************************************************

*************************************************************/

#ifndef LinskerSynapse_H_
#define LinskerSynapse_H_

#ifndef __AVR_ARCH__
#include <type_traits>
#include <concepts>

#include "IntegratorConcept.h"
#include "NeuronConcept.h"
#endif  //__AVR_ARCH__

#include "NormalizableSynapseConcept.h"
#include "SynapseWeightNormalizer.h"
#include "IntegratedSystemWrapper.h"
#include "LinskerSynapseModel.h"
#include "SerializableWrapper.h"
#include "SystemWrapper.h"
#include <cmath>

/**
* Implements a synapse based on (Linsker, 1986)
*/

template <typename TNode1, typename TNode2, typename TIntegrator, typename precission = double>
requires NeuronConcept<TNode1> && NeuronConcept<TNode2> &&
    IntegratorConcept<TIntegrator, SerializableWrapper<
          SystemWrapper<LinskerSynapseModel<precission> > > >

class LinskerSynapse
      : public SerializableWrapper<
            SystemWrapper<LinskerSynapseModel<precission> > > {
 private:
  #ifndef __AVR_ARCH__
    static_assert(std::is_floating_point<precission>::value);
  #endif  //__AVR_ARCH__

  precission m_release_time;

  TNode1 const &m_n1;
  TNode2 &m_n2;

  precission m_last_value_pre;

  const typename TNode1::variable m_n1_variable;
  const typename TNode2::variable m_n2_variable;


  typedef SerializableWrapper<
      SystemWrapper<LinskerSynapseModel<precission> > > System;

  const int m_steps;


 public:
  typedef typename System::precission_t precission_t;
  typedef typename System::variable variable;
  typedef typename System::parameter parameter;
  typedef typename System::ConstructorArgs ConstructorArgs;
  using Normalizer = SynapseWeightNormalizer<TNode2, LinskerSynapse<TNode1, TNode2, TIntegrator, precission>>;

  LinskerSynapse (TNode1 const &n1, typename TNode1::variable v1,
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
          Normalizer::get_instance().add_synapse(this, &m_n2);
        }


  LinskerSynapse (TNode1 const &n1, typename TNode1::variable v1,
                            TNode2 &n2, typename TNode2::variable v2,
                            ConstructorArgs &&args, int steps)
      : m_n1(n1),
        m_n2(n2),
        m_n1_variable(v1),
        m_n2_variable(v2),
        System(args),
        m_steps(steps) {
          for(int i = 0; i < System::n_variables; i++) {
            System::m_variables[i] = 0;
          }
          Normalizer::get_instance().add_synapse(this, &m_n2);
        }

  LinskerSynapse (TNode1 const &n1, TNode2 &n2,
                            LinskerSynapse  const &synapse)
      : m_n1(n1),
        m_n2(n2),
        m_n1_variable(synapse.m_n1_variable),
        m_n2_variable(synapse.m_n2_variable),
        m_steps(synapse.m_steps),
        System(synapse) {
          Normalizer::get_instance().add_synapse(this, &m_n2);
        }

  ~LinskerSynapse() {
      Normalizer::get_instance().remove_synapse(this, &m_n2);
  }

  void step(precission h) {

    System::m_parameters[System::v_pre] = m_n1.get(m_n1_variable);
    precission v_post = m_n2.get(m_n2_variable);
    System::m_parameters[System::v_post] = v_post;
    
    for (int i = 0; i < m_steps; ++i) {
      TIntegrator::step(*this, h, System::m_variables, System::m_parameters);
    }
    Normalizer::get_instance().normalize_weights(this, &m_n2);
    System::m_parameters[System::i] = System::m_variables[System::w] * System::m_parameters[System::v_pre];

  }

  void step(precission h, precission vpre, precission vpost) {

    System::m_parameters[System::v_pre] = vpre;
    precission v_post = vpost;
    System::m_parameters[System::v_post] = v_post;
    
    for (int i = 0; i < m_steps; ++i) {
      TIntegrator::step(*this, h, System::m_variables, System::m_parameters);
    }
    Normalizer::get_instance().normalize_weights(this, &m_n2);
    System::m_parameters[System::i] = System::m_variables[System::w] * System::m_parameters[System::v_pre];

    
  }


  precission get_w_max() const {
    return System::m_parameters[System::w_max];
  }
  
  precission get_weight() const {
    return System::m_variables[System::w];
  }

  void set_weight(precission weight) {
    System::m_variables[System::w] = weight;
  }  
  
 private:

 static_assert(NormalizableSynapseConcept<LinskerSynapse<TNode1, TNode2, TIntegrator, precission>>);
};


#endif /*LINSKER_SYNAPSE_H_*/
