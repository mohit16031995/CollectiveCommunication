#ifndef _ALL_REDUCE_SPARSE_H
#define _ALL_REDUCE_SPARSE_H 

#include "executor.h"

#if defined(_C_BIG)
#include "c_allreduce/c_allreduce_big.h"
#elif defined(_C_SMALL)
#include "c_allreduce/c_allreduce_small.h"
#elif defined(_C_RING)
#include "c_allreduce/c_allreduce_ring.h"
#else
#include "c_allreduce/c_allreduce_recdoubling.h"
#endif

template< class Model, class Params, class Sample, class Loader, class Exec>
class AllReduceSparse : public Executor<Model, Params, Sample, Loader, Exec>
{
public:
    AllReduceSparse() : Executor<Model, Params, Sample, Loader, Exec>()
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

      unsigned dim = this->model_->gradient.size;
      size_t maxbytes = sizeof(unsigned) + dim * sizeof(fp_type);
      struct stream *recvbuf = (struct stream *)malloc(maxbytes);

#if defined(_C_BIG)
      c_allreduce_big<IdxType, fp_type>(this->model_->gradient.values, recvbuf, dim, MPI_COMM_WORLD);
#elif defined(_C_SMALL)
      c_allreduce_small<IdxType, fp_type>(this->model_->gradient.values, recvbuf, dim, MPI_COMM_WORLD);
#elif defined(_C_RING)
      c_allreduce_ring<IdxType, fp_type>(this->model_->gradient.values, recvbuf, dim, MPI_COMM_WORLD);
#else
      c_allreduce_recdoubling<IdxType, fp_type>(this->model_->gradient.values, recvbuf, dim, MPI_COMM_WORLD);
#endif

      if(recvbuf->nofitems == dim) {
        // Dense
        for(unsigned i = 0; i < dim; ++i) {
          this->model_->weights.values[i] += ((fp_type *)recvbuf->items)[i];
        }
      } else {
        // Sparse
        const struct s_item<IdxType, fp_type> *values = (const struct s_item<IdxType, fp_type> *)recvbuf->items;
        for(unsigned i = 0; i < recvbuf->nofitems; ++i) {
          this->model_->weights.values[values[i].idx] += values[i].val;
        }

      }

      free(recvbuf);

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
