#ifndef _ALL_REDUCE_BINOMIAL_H
#define _ALL_REDUCE_BINOMIAL_H 

#include "executor.h"
#include "types/sparse_vector.h"

#define MPI_SPARSE_TAG_INDEX 21

template< class Model, class Params, class Sample, class Loader, class Exec>
class AllReduceBinomialSync : public Executor<Model, Params, Sample, Loader, Exec>
{
private:
    void generate_binomial_tree_direct(int &send_to, vector<int> &receive_from) {
      for(int r = 0; r < ceil(log2(this->world_size_)); ++r) {
        int p2r = pow(2,r);

        int peer = this->rank_ + p2r;
        if(peer < this->world_size_ && this->rank_ < p2r) {
          receive_from.push_back(peer);
        }

        peer = this->rank_ - p2r;
        if(this->rank_ >= p2r && this->rank_ < pow(2, r + 1)) {
          send_to = peer;
        }
      }
    }

    void reduce_sparse(const sparse_vector &my_vector_sparse, sparse_vector &reduced_result_sparse)
    {
      vector<int> receive_from;
      int send_to = -1;

      generate_binomial_tree_direct(send_to, receive_from);

      unsigned int length = receive_from.size();
      vector< MPI_Request > requests(length);
      vector< unsigned int > receive_buffers(length);

      for(unsigned int i = 0; i < length; i++) {
        MPI_Irecv(&receive_buffers[i], 1, MPI_UNSIGNED, receive_from[i], MPI_SPARSE_TAG_INDEX, MPI_COMM_WORLD, &requests[i]);
      }

      int pending_request = length;
      while(pending_request > 0) {
        int index;
        MPI_Status status;
        MPI_Waitany(length, &requests[0], &index, &status); // request should be automatically changed to MPI_REQUEST_NULL by Waitany
        if(index == MPI_UNDEFINED) {
          cout << "Unexpected error!" << endl;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
        sparse_vector rec_sv(receive_buffers[index]);
        MPI_Recv(&rec_sv[0], receive_buffers[index], *mpi_sparse_item_type_, receive_from[index], MPI_SPARSE_TAG_INDEX, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        sum_sparse_intofirst(reduced_result_sparse, rec_sv);
        pending_request--;
      }

      if(send_to >= 0) {
        unsigned int cnt = reduced_result_sparse.size();
        MPI_Send(&cnt, 1, MPI_UNSIGNED, send_to, MPI_SPARSE_TAG_INDEX, MPI_COMM_WORLD);
        MPI_Send(&reduced_result_sparse[0], reduced_result_sparse.size(), *mpi_sparse_item_type_, send_to, MPI_SPARSE_TAG_INDEX,  MPI_COMM_WORLD);
      }
    }
public:
    AllReduceBinomialSync() : Executor<Model, Params, Sample, Loader, Exec>()
    {
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
      return this->rank_;
    }

    void GetModel(hazy::util::Clock &communicate_timer)
    {
      // Always get model after having sent the model updates
    }

    void SendModelUpdate(hazy::util::Clock &communicate_timer)
    {
      communicate_timer.Start();

      sparse_vector my_vector_sparse;
      make_sparse(this->model_->gradient, my_vector_sparse);
      sparse_vector reduced_result_sparse = sparse_vector(my_vector_sparse); // TODO Maybe prevent copying of vector

      reduce_sparse(my_vector_sparse, reduced_result_sparse);

      if(this->rank_ == 0)
      {
        hazy::vector::SVector<fp_type> modelUpdate;
        modelUpdate.size = reduced_result_sparse.size();
        modelUpdate.index = new int[reduced_result_sparse.size()];
        modelUpdate.values = new fp_type[reduced_result_sparse.size()];
        for (unsigned int k = 0; k < reduced_result_sparse.size(); ++k) 
        {
          modelUpdate.index[k] = reduced_result_sparse[k].first;
          modelUpdate.values[k] = reduced_result_sparse[k].second;
        }

        hazy::vector::ScaleAndAdd(this->model_->weights, modelUpdate, 1.0);
        delete[] modelUpdate.values;
        delete[] modelUpdate.index;
      }

      MPI_Bcast(this->model_->weights.values, this->model_->weights.size, MPI_FP_TYPE, 0, MPI_COMM_WORLD);

      communicate_timer.Pause();
    }

    void PreEpoch(hazy::util::Clock &communicate_timer)
    { }

    void PostEpoch(hazy::util::Clock &communicate_timer)
    { }

    void InitStrategy()
    { }

    ~AllReduceBinomialSync()
    {
      MPI_Type_free(mpi_sparse_item_type_);
    }

    MPI_Datatype *mpi_sparse_item_type_;
};

#endif
