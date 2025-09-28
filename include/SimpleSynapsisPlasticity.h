/*************************************************************


*************************************************************/

#ifndef SIMPLE_SYNAPSIS_PLASTICITY_H_
#define SIMPLE_SYNAPSIS_PLASTICITY_H_

#ifndef __AVR_ARCH__
#include <type_traits>

#include "IntegratorConcept.h"
#include "NeuronConcept.h"
#endif  //__AVR_ARCH__

#include "SimpleSynapsisPlasticityModel.h"
#include "IntegratedSystemWrapper.h"
#include "SerializableWrapper.h"
#include "SystemWrapper.h"
#include <cmath>

/**
 * @brief Implements a simple synapsis with plasticity
 */
template <typename TNode1, typename TNode2, typename TIntegrator,
          typename precission = double>
requires NeuronConcept<TNode1> && NeuronConcept<TNode2> &&
    IntegratorConcept<TIntegrator, SerializableWrapper<
          SystemWrapper<SimpleSynapsisPlasticityModel<precission> > > >

class SimpleSynapsisPlasticity
    : public SerializableWrapper<
          SystemWrapper<SimpleSynapsisPlasticityModel<precission> > >{

 private:
#ifndef __AVR_ARCH__
  static_assert(std::is_floating_point<precission>::value);
#endif // __AVR_ARCH__

  precission m_release_time;

  TNode1 const &m_n1;
  TNode2 &m_n2;

  precission m_last_value_pre;

  const typename TNode1::variable m_n1_variable;
  const typename TNode2::variable m_n2_variable;


  typedef SerializableWrapper<
      SystemWrapper<SimpleSynapsisPlasticityModel<precission> > > System;

  const int m_steps;

 public:
  typedef typename System::precission_t precission_t;
  typedef typename System::variable variable;
  typedef typename System::parameter parameter;
  typedef typename System::ConstructorArgs ConstructorArgs;

  SimpleSynapsisPlasticity(TNode1 const &n1, typename TNode1::variable v1, TNode2 &n2, typename TNode2::variable v2, ConstructorArgs &args, int steps)
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


  SimpleSynapsisPlasticity(TNode1 const &n1, typename TNode1::variable v1,
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
        }

  SimpleSynapsisPlasticity(TNode1 const &n1, TNode2 &n2,
                            SimpleSynapsisPlasticity const &synapse)
      : m_n1(n1),
        m_n2(n2),
        m_n1_variable(synapse.m_n1_variable),
        m_n2_variable(synapse.m_n2_variable),
        m_steps(synapse.m_steps),
        System(synapse) {}


  void step(precission h) {

    System::m_parameters[System::v_pre] = m_n1.get(m_n1_variable);
    precission v_post = m_n2.get(m_n2_variable);

    for (int i = 0; i < m_steps; ++i) {
      TIntegrator::step(*this, h, System::m_variables, System::m_parameters);
    }

    /* Degradation calculation begin*/
    if (System::m_parameters[System::v_pre] >= System::m_parameters[System::Vfast]) { 
    /* Synapsis is active*/
    System::m_parameters[System::steps_remaining] = System::m_parameters[System::steps_for_degradation];
    
      if (System::m_parameters[System::previous_System::m_parameters[System::v_pre]] < System::m_parameters[System::Vfast]) {
          System::m_parameters[System::gfast] += 0.1; // LTP (Long-term Potentiation)
          if (System::m_parameters[System::gfast] > 1) {
              System::m_parameters[System::gfast] = 1; // max gfast
          }
      }
    } else { 
      /*Synapsis is not active*/
      if (System::m_parameters[System::steps_remaining] == 0) {
          System::m_parameters[System::gfast] -= 0.1; // LTD (Long-term Depression)
          if (System::m_parameters[System::gfast] < 0) {
              System::m_parameters[System::gfast] = 0; // min gfast
          }
      } else {
          System::m_parameters[System::steps_remaining] -= 1;
      }
    }
    System::m_parameters[System::previous_v_pre] = System::m_parameters[System::v_pre];
    /* Degradation calculation end*/



    System::m_parameters[System::ifast] = (System::m_parameters[System::gfast] * (v_post - System::m_parameters[System::Esyn])) /
                                          (1 + exp(System::m_parameters[System::sfast] * (System::m_parameters[System::Vfast] - System::m_parameters[System::v_pre])));

    /* islow will be ignored*/

    System::m_parameters[System::i] = System::m_parameters[System::ifast];

  }

  void step(precission h, precission vpre, precission vpost) {

    System::m_parameters[System::v_pre] = vpre;
    precission v_post = vpost;

    for (int i = 0; i < m_steps; ++i) {
      TIntegrator::step(*this, h, System::m_variables, System::m_parameters);
    }

    /* Degradation calculation begin*/
    if (vpre >= System::m_parameters[System::Vfast]) { 
      /* Synapsis is active*/
      System::m_parameters[System::steps_remaining] = System::m_parameters[System::steps_for_degradation];
    
      if (System::m_parameters[System::previous_vpre] != -999999 && System::m_parameters[System::previous_vpre] < System::m_parameters[System::Vfast]) {
          System::m_parameters[System::gfast] += 0.20; // LTP (Long-term Potentiation)
          if (System::m_parameters[System::gfast] > 1) {
              System::m_parameters[System::gfast] = 1; // max gfast
          }
      }
    } else { 
      /*Synapsis is not active*/
      if (System::m_parameters[System::steps_remaining] == 0) {
          System::m_parameters[System::steps_remaining] = System::m_parameters[System::steps_for_degradation];
          System::m_parameters[System::gfast] -= 0.005; // LTD (Long-term Depression)
          if (System::m_parameters[System::gfast] < 0) {
              System::m_parameters[System::gfast] = 0; // min gfast
          }
      } else {
          System::m_parameters[System::steps_remaining] -= 1;
      }
    }
    System::m_parameters[System::previous_vpre] = vpre;
    /* Degradation calculation end*/

    System::m_parameters[System::ifast] = (System::m_parameters[System::gfast] * (v_post - System::m_parameters[System::Esyn])) /
                                          (1 + exp(System::m_parameters[System::sfast] * (System::m_parameters[System::Vfast] - System::m_parameters[System::v_pre])));

    /* islow will be ignored*/

    System::m_parameters[System::i] = System::m_parameters[System::ifast];

  }
};

#endif /*SIMPLE_SYNAPSIS_PLASTICITY_H_*/
