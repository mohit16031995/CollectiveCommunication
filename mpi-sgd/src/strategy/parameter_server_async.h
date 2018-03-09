#ifndef _PARAMETER_SERVER_ASYNC_H
#define _PARAMETER_SERVER_ASYNC_H

#include "executor.h"

#define MPI_TAG_GET_MODEL 11
#define MPI_TAG_SEND_UPDATE 12
#define MPI_TAG_NOTIFY_END_EPOCH 13

template< class Model, class Params, class Sample, class Loader, class Exec>
class ParameterServerAsync : public Executor<Model, Params, Sample, Loader, Exec>
{
public:
    ParameterServerAsync() : Executor<Model, Params, Sample, Loader, Exec>()
    {
      this->runThreadPool_ = this->rank_ != 0;
    }

    int GetWorkerNumber()
    {
      if(this->rank_ == 0)
        return -1;

      return this->rank_-1;
    }

    void InitStrategy()
    {
      if(this->rank_ == 0)
      {
        // Prepare requests and stati array
        this->req_ = new MPI_Request[(this->world_size_ -1)*2];
        this->stat_ = new MPI_Status[(this->world_size_ -1)*2];

        this->modelUpdates_ = new fp_type*[this->world_size_ -1];
        
        // Set receive and issend
        for(int i = 0; i < this->world_size_ -1; ++i)
        {
          this->modelUpdates_[i] = new fp_type[this->model_->weights.size];
          MPI_Irecv(modelUpdates_[i], this->model_->weights.size, MPI_FP_TYPE, i+1, MPI_ANY_TAG, MPI_COMM_WORLD, &req_[i*2]);
          MPI_Issend(this->model_->weights.values, this->model_->weights.size, MPI_FP_TYPE, i+1, MPI_TAG_GET_MODEL, MPI_COMM_WORLD, &req_[(i*2)+1]);
        }
      }
    }

    void GetModel(hazy::util::Clock &communicate_timer)
    {
      if(this->rank_ != 0)
      {
        communicate_timer.Start();
        MPI_Status stat;
        MPI_Recv(this->model_->weights.values, this->model_->weights.size, MPI_FP_TYPE, 0, MPI_TAG_GET_MODEL, MPI_COMM_WORLD, &stat);
        communicate_timer.Pause();
      }
    }

    void SendModelUpdate(hazy::util::Clock &communicate_timer)
    {
      if(this->rank_ != 0)
      {
        communicate_timer.Start();
        MPI_Send(this->model_->gradient.values, this->model_->gradient.size, MPI_FP_TYPE, 0, MPI_TAG_SEND_UPDATE , MPI_COMM_WORLD);
        communicate_timer.Pause();
      }
    }

    void PreEpoch(hazy::util::Clock &communicate_timer)
    {
      if(this->rank_ == 0)
      {
        communicate_timer.value = 0.0;
      }
    }

    void PostEpoch(hazy::util::Clock &communicate_timer)
    {
      if(this->rank_ != 0)
      {
        communicate_timer.Start();
        // Notify server
        MPI_Send(NULL, 0, MPI_BYTE, 0, MPI_TAG_NOTIFY_END_EPOCH, MPI_COMM_WORLD);
        communicate_timer.Pause();
      }
      else
      {
        // Handle all request until everyone has notified a post epoch
        int workerFinished = 0;
        while(workerFinished < this->world_size_ -1)
        {
          int count = 0;
          int index[(this->world_size_-1)*2];
          MPI_Waitsome((this->world_size_-1)*2, req_, &count, index, stat_);

          for(int i = 0; i < count; ++i)
          {
            int sender = 0;
            if(index[i] % 2 == 0)
            {
              // Model update was received
              sender = index[i]/2;

              switch(stat_[i].MPI_TAG)
              {
                case MPI_TAG_SEND_UPDATE :
                  {
                    hazy::vector::FVector<fp_type> modelUpdate(modelUpdates_[sender], this->model_->weights.size);

                    hazy::vector::ScaleAndAdd(this->model_->weights, modelUpdate, 1.0);
                  }
                  break;
                case MPI_TAG_NOTIFY_END_EPOCH :
                  workerFinished++;
                  break;
                default :
                  printf("Received message with unknown tag\n");
                  MPI_Abort(MPI_COMM_WORLD, 0);
              }

              MPI_Irecv(modelUpdates_[sender], this->model_->weights.size, MPI_FP_TYPE, sender+1, MPI_ANY_TAG, MPI_COMM_WORLD, &req_[index[i]]);
            }
            else
            {
              sender = (index[i]-1)/2;

              MPI_Issend(this->model_->weights.values, this->model_->weights.size, MPI_FP_TYPE, sender+1, MPI_TAG_GET_MODEL, MPI_COMM_WORLD, &req_[index[i]]);
            }
          }
        }
      }

      if(this->rank_ != 0)
        communicate_timer.Start();

      MPI_Barrier(MPI_COMM_WORLD);

      if(this->rank_ != 0)
        communicate_timer.Pause();
    }

    ~ParameterServerAsync()
    {
      if(this->rank_ == 0)
      {
        // Cleanup statis and requests
        delete[] this->req_;
        delete[] this->stat_;

        // Cleanup modelUpdate buffers
        for(int i = 0; i < this->world_size_ - 1; ++i)
        {
          delete[] modelUpdates_[i];
        }
        delete[] modelUpdates_;
      }
    }

    MPI_Request *req_;
    MPI_Status *stat_;
    fp_type** modelUpdates_;
};

#endif
