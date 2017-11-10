#include "typedefs.h"

using namespace std;

typedef map<Dart_handle, Dart_handle> bijection;

template<class LCC>
class LccComparer
{
private:
  LCC lcc1, lcc2;
  bijection mF; // The traversal between the two LCC's

public:
  bool is_isomorphic_function(LCC *lcc1, LCC *lcc2, bijection &f)
  {
    for(typename LCC::Dart_range::iterator it = lcc1->darts().begin(); it != lcc1->darts().end(); ++it )
    {
      for (int i = 0; i <= lcc1->dimension; i++)
      {
        if (f[lcc1->beta(it, i)] != lcc2->beta(f[it], i))
          return false;
      }
    }

    return true;
  }

  void traverseAndBuildMatching(LCC *lcc1, LCC *lcc2, Dart_handle d1, Dart_handle d2)
  {
    mF.clear();

    for(typename LCC::Dart_range::iterator it = lcc1->darts().begin(); it != lcc1->darts().end(); ++it )
    {
      mF[it] = NULL;
    }

    mF[d1] = d2;

    vector<Dart_handle> s;
    s.push_back(d1);

    while (s.size() > 0)
    {
      Dart_handle d = s.back();
      s.pop_back();
      for (int i = 0; i <= lcc1->dimension; i++)
      {
        if (lcc1->beta(d, i) != lcc1->null_handle && mF[lcc1->beta(d, i)] == NULL)
        {
          mF[lcc1->beta(d, i)] = lcc2->beta(mF[d], i);
          s.push_back(lcc1->beta(d, i));
        }
      }
    }

    mF[lcc1->null_handle] = lcc2->null_handle;
  }

  bool compare(LCC &lcc1, LCC &lcc2)
  {
    Dart_handle d0 = lcc1.darts().begin();
    for(typename LCC::Dart_range::iterator it = lcc2.darts().begin(); it != lcc2.darts().end(); ++it )
    {
      traverseAndBuildMatching(&lcc1, &lcc2, d0, it);
      if (is_isomorphic_function(&lcc1, &lcc2, mF))
        return true;
    }

    return false;
  }
};
