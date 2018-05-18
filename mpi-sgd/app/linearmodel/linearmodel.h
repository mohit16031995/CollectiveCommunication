#ifndef _LINEARMODEL_H
#define _LINEARMODEL_H

#include "hazy/vector/fvector.h"
#include "hazy/vector/svector.h"
#include "hazy/vector/scale_add-inl.h"
#include "hazy/vector/operations-inl.h"

struct LinearModel
{
  hazy::vector::FVector<fp_type> weights;
  hazy::vector::FVector<fp_type> gradient;
  hazy::vector::FVector<fp_type> * local_gradients;

  float batch_step_size;    //!< current batch step size
  DECREASING_STEPSIZES_ONLY(unsigned long long k;); //!< number of processed element for batch_step_size calculation

  size_t nthreads_;

  explicit LinearModel(unsigned dim, unsigned nthreads)
  {
    nthreads_ = nthreads;
    weights.values = new fp_type[dim];
    gradient.values = new fp_type[dim];
    weights.size = dim;
    gradient.size = dim;

    for(unsigned i = 0; i < dim; ++i)
    {
      weights.values[i] = 0;
      gradient.values[i] = 0;
    }

    if(nthreads > 1)
    {
      local_gradients = new hazy::vector::FVector<fp_type>[nthreads];
    }
    else
    {
      local_gradients = &gradient;
    }
  }

  void initLocalVars(unsigned dim, unsigned tid)
  {
    local_gradients[tid].values = (fp_type*)aligned_alloc(CACHE_LINE_SIZE, dim * sizeof(fp_type));
    local_gradients[tid].size = dim;
    for (unsigned i = 0; i < dim; ++i) {
      local_gradients[tid].values[i] = 0;
    }
  }

  ~LinearModel()
  {
    if(nthreads_ > 1)
    {
      for(size_t j = 0; j < nthreads_; j++) {
        free(local_gradients[j].values);
      }
      delete[] local_gradients;
    }
    delete[] weights.values;
    delete[] gradient.values;
  }
};

struct LinearModelParams
{
  unsigned batch_size;        //!< batch size
  float step_size;            //!< stepsize (decayed by step decay at each epoch and/or after each mini batch...)
  float step_decay;           //!< factor to modify step_size by each epoch
  float beta;                 //!< the exponent constant for the stepsizes
  float lasso_regularizer;    //!< lasso regularizer hyperparamter
  uint64_t numSamplesProc;    //!< Number of samples on current worker
  uint64_t totalNumSamples;   //!< Total number of samples
  uint64_t maxSamplesProc;    //!< Max number on any worker
  unsigned ndim;              //!< number of features, length of degrees
  unsigned nthreads;          //!< Number of threads on working processes
  unsigned quantizationLevel; //!< Quantization level (only affects LinReg and if quantization is > 0)
  unsigned quantization;      //!< 0: No quantization / 1: Quatize samples / 2: quantize gradients / 3: quantize samples and gradient / 4: Quantize all Model / Sample and Gradient)
  hazy::vector::FVector<fp_type> *x_hat;

  LinearModelParams(fp_type stepsize, fp_type stepdecay) : step_size(stepsize), step_decay(stepdecay) { }
};

struct LinearModelSample
{
  fp_type value;            //!< rating of this example
  SPARSE_ONLY(hazy::vector::SVector<const fp_type> vector;); //!< feature vector
  DENSE_ONLY( hazy::vector::FVector<const fp_type> vector;); //!< feature vector

  LinearModelSample()
  { }

#ifdef _SPARSE
  LinearModelSample(
    fp_type val,
    fp_type const *values,
    int *index, 
    unsigned len,
    int dimension
  ) : value(val), vector(values, index, len) { }
#else
  LinearModelSample(
    fp_type val,
    fp_type const *values,
    int *index, 
    unsigned len,
    int dimension
  ) : value(val) {

    hazy::vector::SVector<const fp_type> sparse_vector(values, index, len);

    fp_type *zeros = new fp_type[dimension];
    hazy::vector::FVector<fp_type> *temp_vector = new hazy::vector::FVector<fp_type>(zeros, dimension);

    // set to zero
    hazy::vector::Zero(*temp_vector);
    
    hazy::vector::ScaleAndAdd(
		      *temp_vector,
		      sparse_vector,
		      1.0
    );

    vector.size = (*temp_vector).size;
    vector.values = (*temp_vector).values;
  }
#endif

  LinearModelSample(const LinearModelSample &o) {
    value = o.value;
    vector.values = o.vector.values;
    SPARSE_ONLY(vector.index = o.vector.index;);
    vector.size = o.vector.size;
  }

};

#endif
