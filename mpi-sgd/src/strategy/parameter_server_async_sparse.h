#ifndef _PARAMETER_SERVER_ASYNC_SPARSE_H
#define _PARAMETER_SERVER_ASYNC_SPARSE_H

#include "executor.h"
#include "types/sparse_vector.h"

#define MPI_TAG_GET_MODEL 11
#define MPI_TAG_SEND_UPDATE 12
#define MPI_TAG_NOTIFY_END_EPOCH 13

template< class Model, class Params, class Sample, class Loader, class Exec>
class ParameterServerAsyncSparse : public Executor<Model, Params, Sample, Loader, Exec>
{
public:
    ParameterServerAsyncSparse() : Executor<Model, Params, Sample, Loader, Exec>()
    {
      this->runThreadPool_ = this->rank_ != 0;

      mpi_sparse_item_type_ = new MPI_Datatype;
      const int nitems = 2; // Pair of index / data_type
      int blocklengths[nitems] = {1, 1};
      MPI_Datatype types[nitems] = {MPI_UNSIGNED, MPI_FP_TYPE};

      MPI_Aint offsets[nitems];
      offsets[0] = offsetof(sparse_item, first);
      offsets[1] = offsetof(sparse_item, second);

      MPI_Type_create_struct(nitems, blocklengths, offsets, types, mpi_sparse_item_type_);
      MPI_Type_commit(mpi_sparse_item_type_);
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

        this->modelUpdatesSize_.resize(this->world_size_-1);
        // Set receive and issend
        for(int i = 0; i < this->world_size_ -1; ++i)
        {
          MPI_Irecv(&modelUpdatesSize_[i], 1, MPI_UNSIGNED, i+1, MPI_ANY_TAG, MPI_COMM_WORLD, &req_[i*2]);
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
        sparse_vector v;
        make_sparse(this->model_->gradient, v);
        unsigned int cnt = v.size();
        MPI_Send(&cnt, 1, MPI_UNSIGNED, 0, MPI_TAG_SEND_UPDATE , MPI_COMM_WORLD);
        MPI_Send(&v[0], v.size(), *this->mpi_sparse_item_type_, 0, MPI_TAG_SEND_UPDATE , MPI_COMM_WORLD);
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
                    sparse_vector rec_sv(modelUpdatesSize_[sender]);
                    MPI_Recv(&rec_sv[0], modelUpdatesSize_[sender], *mpi_sparse_item_type_, sender+1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    hazy::vector::SVector<fp_type> modelUpdate;
                    modelUpdate.size = rec_sv.size();
                    modelUpdate.index = new int[rec_sv.size()];
                    modelUpdate.values = new fp_type[rec_sv.size()];
                    for (unsigned int k = 0; k < rec_sv.size(); ++k) 
                    {
                      modelUpdate.index[k] = rec_sv[k].first;
                      modelUpdate.values[k] = rec_sv[k].second;
                    }

                    hazy::vector::ScaleAndAdd(this->model_->weights, modelUpdate, 1.0);
                    delete[] modelUpdate.values;
                    delete[] modelUpdate.index;
                  }
                  break;
                case MPI_TAG_NOTIFY_END_EPOCH :
                  workerFinished++;
                  break;
                default :
                  printf("Received message with unknown tag\n");
                  MPI_Abort(MPI_COMM_WORLD, 0);
              }

              MPI_Irecv(&modelUpdatesSize_[sender], 1, MPI_UNSIGNED, sender+1, MPI_ANY_TAG, MPI_COMM_WORLD, &req_[index[i]]);
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

    ~ParameterServerAsyncSparse()
    {
      if(this->rank_ == 0)
      {
        // Cleanup statis and requests
        delete[] this->req_;
        delete[] this->stat_;
      }
      MPI_Type_free(mpi_sparse_item_type_);
    }

    MPI_Request *req_;
    MPI_Status *stat_;
    MPI_Datatype *mpi_sparse_item_type_;
    vector< unsigned int > modelUpdatesSize_;
};

#endif
