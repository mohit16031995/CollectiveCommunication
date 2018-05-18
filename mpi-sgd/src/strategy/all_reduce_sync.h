#ifndef _ALL_REDUCE_SYNC_H
#define _ALL_REDUCE_SYNC_H 

#include "executor.h"

template< class Model, class Params, class Sample, class Loader, class Exec>
class AllReduceSync : public Executor<Model, Params, Sample, Loader, Exec>
{
public:
    AllReduceSync() : Executor<Model, Params, Sample, Loader, Exec>()
    { }

    int GetWorkerNumber()
    {
      return this->rank_;
    }

    void GetModel(hazy::util::Clock &communicate_timer)
    {
      // Always get model after having sent the model updates
    }

    void SendModelUpdate(hazy::util::Clock &communicate_timer)
    {
      communicate_timer.Start();

      fp_type *d = new fp_type[this->model_->gradient.size];
      hazy::vector::FVector<fp_type> gradients(d, this->model_->gradient.size);
      hazy::vector::Zero(gradients);
      MPI_Allreduce(this->model_->gradient.values, gradients.values, this->model_->gradient.size, MPI_FP_TYPE, MPI_SUM, MPI_COMM_WORLD);
      hazy::vector::ScaleAndAdd(this->model_->weights, gradients, 1.0);
      delete[] d;

      communicate_timer.Pause();
    }

    void PreEpoch(hazy::util::Clock &communicate_timer)
    { }

    void PostEpoch(hazy::util::Clock &communicate_timer)
    { }

    void InitStrategy()
    { }
};

#endif
