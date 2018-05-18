#ifndef _PARAMETER_SERVER_SYNC_H
#define _PARAMETER_SERVER_SYNC_H

#include "executor.h"

template< class Model, class Params, class Sample, class Loader, class Exec>
class ParameterServerSync : public Executor<Model, Params, Sample, Loader, Exec>
{
public:
    ParameterServerSync() : Executor<Model, Params, Sample, Loader, Exec>()
    {
      this->runThreadPool_ = this->rank_ != 0;
    }

    int GetWorkerNumber()
    {
      if(this->rank_ == 0)
        return -1;

      return this->rank_-1;
    }

    void GetModel(hazy::util::Clock &communicate_timer)
    {
      if(this->rank_ != 0)
        communicate_timer.Start();

      MPI_Bcast(this->model_->weights.values, this->model_->weights.size, MPI_FP_TYPE, 0, MPI_COMM_WORLD);

      if(this->rank_ != 0)
        communicate_timer.Pause();
    }

    void SendModelUpdate(hazy::util::Clock &communicate_timer)
    {
      if(this->rank_ != 0)
        communicate_timer.Start();

      fp_type *d = new fp_type[this->model_->gradient.size];
      hazy::vector::FVector<fp_type> gradients(d, this->model_->gradient.size);
      hazy::vector::Zero(gradients);
      MPI_Reduce(this->model_->gradient.values, gradients.values, this->model_->gradient.size, MPI_FP_TYPE, MPI_SUM, 0, MPI_COMM_WORLD);

      if(this->rank_ == 0)
      {
        hazy::vector::ScaleAndAdd(this->model_->weights, gradients, 1.0);
      }

      delete[] d;

      if(this->rank_ != 0)
        communicate_timer.Pause();
    }

    void PreEpoch(hazy::util::Clock &communicate_timer)
    {
      if(this->rank_ == 0)
        communicate_timer.value = 0.0;
    }

    void PostEpoch(hazy::util::Clock &communicate_timer)
    { }

    void InitStrategy()
    { }
};

#endif
