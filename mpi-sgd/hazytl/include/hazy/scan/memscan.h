// Copyright 2012 Chris Re, Victor Bittorf
//
 //Licensed under the Apache License, Version 2.0 (the "License");
 //you may not use this file except in compliance with the License.
 //You may obtain a copy of the License at
 //    http://www.apache.org/licenses/LICENSE-2.0
 //Unless required by applicable law or agreed to in writing, software
 //distributed under the License is distributed on an "AS IS" BASIS,
 //WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 //See the License for the specific language governing permissions and
 //limitations under the License.

// The Hazy Project, http://research.cs.wisc.edu/hazy/
// Author : Victor Bittorf (bittorf [at] cs.wisc.edu)

#ifndef HAZY_HOGWILD_FILE_SCAN_H
#define HAZY_HOGWILD_FILE_SCAN_H

#include "hazy/util/simple_random-inl.h"

#include "hazy/vector/fvector.h"
#include "hazy/scan/sampleblock.h"

namespace hazy {
namespace scan {

/*! \brief A simple scanner that permutes examples stored in memroy.
 * Returns all examples in a single page
 */
template <class Sample>
class MemoryScan {
 public:
  /*! \brief Makes a new scanner over the given vector of examples
   */
  MemoryScan(vector::FVector<Sample> &fv) { 
    blk_.ex.size = fv.size;
    blk_.ex.values = fv.values;
    blk_.perm.size = 0;
    blk_.perm.values = NULL;
    has_next_ = true;
  }

  ~MemoryScan() {
    if (blk_.perm.values != NULL) {
      delete [] blk_.perm.values;
    }
  }
  
  /*! \brief returns true if there it is valid to call Next()
   */
  bool HasNext() { return has_next_; }

  /*! \brief Gets the next block of examples
   * First permutes the block of examples and the returns the block
   */
  SampleBlock<Sample>& Next() {

    size_t size = blk_.ex.size;

    if (blk_.perm.values != NULL) {
      delete [] blk_.perm.values;
      blk_.perm.values = NULL;
    }
    blk_.perm.size = size;
    blk_.perm.values = new size_t[size];
    for (size_t i = 0; i < size; i++) {
      blk_.perm.values[i] = i;
    }

    util::SimpleRandom &rand = util::SimpleRandom::GetInstance();
    rand.LazyPODShuffle(blk_.perm.values, size);
    has_next_ = false;
    return blk_;
  }

  /*! \brief Resets the scanner to the begining.
   */
  void Reset() { has_next_ = true; }

 private:
  SampleBlock<Sample> blk_;
  bool has_next_;
};


template <class Sample>
class MemoryScanNoPermutation {
 public:
  /*! \brief Makes a new scanner over the given vector of examples
   */
  MemoryScanNoPermutation(vector::FVector<Sample> &fv) { 
    blk_.ex.size = fv.size;
    blk_.ex.values = fv.values;
    blk_.perm.size = 0;
    blk_.perm.values = NULL;
    has_next_ = true;
  }

  ~MemoryScanNoPermutation() {
    if (blk_.perm.values != NULL) {
      delete [] blk_.perm.values;
    }
  }
  
  /*! \brief returns true if there it is valid to call Next()
   */
  bool HasNext() { return has_next_; }

  /*! \brief Gets the next block of examples
   * Just returns the block
   */
  SampleBlock<Sample>& Next() {

    size_t size = blk_.ex.size;

    if (blk_.perm.values != NULL) {
      delete [] blk_.perm.values;
      blk_.perm.values = NULL;
    }
    blk_.perm.size = size;
    blk_.perm.values = new size_t[size];
    for (size_t i = 0; i < size; i++) {
      blk_.perm.values[i] = i;
    }

    has_next_ = false;
    return blk_;
  }

  /*! \brief Resets the scanner to the begining.
   */
  void Reset() { has_next_ = true; }

 private:
  SampleBlock<Sample> blk_;
  bool has_next_;
};

} // namespace hogwild
} // namespace hazy
#endif
