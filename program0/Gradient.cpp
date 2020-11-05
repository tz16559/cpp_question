

using namespace std;
double raw_grad_finder(double x1[], double y1[], double x2[], double y2[], uint k, uint kp){
  // find the gradient of a straight line
  // here usnig this function to find the gradients among different combinations of points forming straight lines
  return (y2[kp]-y1[k])/(x2[kp]-x1[k]);
}

vector<pair<double,uint> > array_sorter(double arr[], uint n){
  // sort an array, together with its changes to the indices of the elements after sorting
  vector<pair<double, uint> > vp;
  for (uint i=0; i<n; i++) {vp.push_back(make_pair(arr[i], i));}
  sort(vp.begin(), vp.end());
  return vp;
}

uint * closest_finder(double arr1[], double arr2[], double arr3[], double arr4[], uint n){
  // find the only element in an array that has the closest value to the value of an element in each of the other arrays
  // here the function is designed for four arrays
  // using this function to find the closest values in gradients from different conbinations for each coordinate set k
  double diff=DBL_MAX;
  uint m0=0, m1=0, m2=0, m3=0;
  static uint min_idx[] = {0,0,0,0};

  while (m0<n && m1<n && m2<n && m3<n){
    double minimum = min(arr1[m0], min(arr2[m1], min(arr3[m2],arr4[m3])));
    double maximum = max(arr1[m0], max(arr2[m1], max(arr3[m2],arr4[m3])));

      if (maximum-minimum < diff){
        min_idx[0] = m0, min_idx[1] = m1, min_idx[2] = m2, min_idx[3] = m3;
        diff = maximum-minimum;
      }

      if (diff == 0) break;

      if      (arr1[m0] == minimum) m0++;
      else if (arr2[m1] == minimum) m1++;
      else if (arr3[m2] == minimum) m2++;
      else                       m3++;
  }
  return min_idx;
}
