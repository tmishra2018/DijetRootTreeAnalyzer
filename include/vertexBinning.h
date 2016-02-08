#pragma once

#include <cmath>
#include <vector>
#include <utility>

class VertexBinning {
  public:
    VertexBinning() {
      fillVertexBins();
    }

    int getVertexBin(int n) {
      std::vector<std::pair<int, int> >::const_iterator it = mVertexBins.begin();
      for (; it != mVertexBins.end(); ++it) {
        std::pair<int, int> bin = *it;
        if (n >= bin.first && n < bin.second) {
          return it - mVertexBins.begin();
        }
      }

      return -1;
    }

    size_t size() const {
      return mVertexBins.size();
    }

    std::pair<int, int> getBinValue(int bin) const {
      return mVertexBins[bin];
    }

    std::vector<std::pair<int, int> > getBinning(int n = -1) const {
      if (n < 0) {
        n = size();
      }
      return std::vector<std::pair<int, int> >(mVertexBins.begin(), mVertexBins.begin() + n);
    }

    std::vector<std::pair<int, int> > getBinning(unsigned int from, unsigned int to) const {
      if (to > size()) {
        to = size();
      }

      return std::vector<std::pair<int, int> >(mVertexBins.begin() + from, mVertexBins.begin() + to);
    }

  private:
    std::vector<std::pair<int, int> > mVertexBins;
    int step = 1 ;
    int nVtx_max= 50;

    void fillVertexBins() {

      int n_bin = nVtx_max/ step;
      
      for(int jj=0; jj<n_bin ; jj++){	 
	int vtx_min = ( (jj-1)*step)+step;
	int vtx_max = (jj*step)+step;
	
	mVertexBins.push_back(std::make_pair(vtx_min, vtx_max));
	
      }
            
      //federico
      //      mVertexBins.push_back(std::make_pair(0, 6)); // 6 escluso : Npv minore = 5
      //      mVertexBins.push_back(std::make_pair(6, 13));
      //      mVertexBins.push_back(std::make_pair(13, 35));

    }

};
