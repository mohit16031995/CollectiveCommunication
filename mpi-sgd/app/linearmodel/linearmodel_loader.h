#ifndef _LINEARMODEL_LOADER_H
#define _LINEARMODEL_LOADER_H

#include "hazy/vector/dot-inl.h"
#include "hazy/vector/operations-inl.h"
#include "hazy/types/entry.h"

class LinearModelLoader
{
public:
  template <class Scan> static size_t LoadSamples(Scan &scan, hazy::vector::FVector<LinearModelSample> &ex, int dimension)
  {
    std::vector<LinearModelSample> examps;
    int lastrow = -1;
    double rating = 0.0;
    std::vector<fp_type> data;
    std::vector<int> index;

    int max_col = 0;

    while (scan.HasNext()) {
      const hazy::types::Entry &e = scan.Next();
      if (lastrow == -1) {
        // this will be the case at the beginning
        lastrow = e.row;
      }
      if ((lastrow != e.row) || (!scan.HasNext())) {
        // finish off the previous vector and start a new one
        lastrow = e.row;

        if(!scan.HasNext()) {
          data.push_back(e.rating);
          index.push_back(e.col);
        }

        fp_type *d = new fp_type[data.size()];
        int *i = new int[data.size()];
        for (size_t j = 0; j < data.size(); j++) {
          d[j] = data[j];
          i[j] = index[j];
        }
        LinearModelSample temp(rating, d, i, data.size(), dimension);
        examps.push_back(temp);
        rating = 0.0;
        data.clear();
        index.clear();
      }

      if (e.col < 0) {
        rating = e.rating;
      } else {
        if (e.col > max_col) {
          max_col = e.col;
        }
        data.push_back(e.rating);
        index.push_back(e.col);
      }
    }

    // Copy from temp vector into persistent memory
    ex.size = examps.size();
    ex.values = new LinearModelSample[ex.size];
    for (size_t i = 0; i < ex.size; i++) {
      new (&ex.values[i]) LinearModelSample(examps[i]);
    }
    return max_col+1;
  }
};

#endif
