#ifndef TRACK_H
#define TRACK_H


#include <vector>
class Track
{
private:
uint p;
// double r[8];
// double t[8];
// double A[p];
// double B[p];
// double C[p];
// double D[p];
// double x_track[p];
// double x_track_err[p];
// double y_track[p];
// double y_track_err[p];
double x_track_p[p];
double y_track_p[p];
double x_track_m[p];
double y_track_m[p];


// double arr[];
// double diff=DBL_MAX;
// double minimum;
// double maximum;
// uint add;
// uint k;
// uint kp;
// uint n;
// uint m0=0;
// uint m1=0;
// uint m2=0;
// uint m3=0;
// static uint min_idx[] = {0,0,0,0};
// vector<pair<double,uint> >  vp;


public:
  Track(double, double, double, uint, uint, uint);
  Track::Track(double t0[],x0[],y[0],uint n0, uint m0, uint add0 ){
    t = t0;
    x = x0;
    y = y0;
    n = n0;
    m = m0;
    add = add0;
  }



  void track_coord(double r0[], double t0[], uint add0);

  Track::track_coord(double r0[], double t0[], uint add0){

    for (uint k=0; k<n; k++){
      x_track_p[k] = track_coord_finder(x, t, k, add0, 1);
      y_track_p[k] = track_coord_finder(y, t, k, add0, 1);
      if (t[k] != t[k+add0]){
        x_track_m[k] = track_coord_finder(x, t, k, add0, -1);
        y_track_m[k] = track_coord_finder(y, t, k, add0, -1);
      }
      else{
        x_track_m[k] = x_track_p[k];
        y_track_m[k] = y_track_p[k];
      }
    }
  }


//   void set_track_coord(double, double, uint, uint, uint);
//   void set_grad_finder(double, double, double, double, uint, uint);
//   double track_coord_finder(){
//     return r[k+add]+c*r[k]*t[k+add]/t[k])/(1+c*t[k+add]/t[k];
//   }
//
//
//   double track_coord_error(){
//     return sqrt((c*(r[k]/t[k])/(1+c*t[k+add]/t[k])-c*1/t[k]*(r[k+add]+c*r[k]*t[k+add]/t[k])/((1+c*t[k+add]/t[k])*(1+c*t[k+add]/t[k])))*
//                 (c*(r[k]/t[k])/(1+c*t[k+add]/t[k])-c*1/t[k]*(r[k+add]+c*r[k]*t[k+add]/t[k])/((1+c*t[k+add]/t[k])*(1+c*t[k+add]/t[k])))*0.25*0.25
//               +(-c*(r[k]*t[k+add]/(t[k]*t[k]))/(1+c*t[k+add]/t[k])+c*t[k+add]/(t[k]*t[k])*(r[k+add]+c*r[k]*t[k+add]/t[k])/((1+c*t[k+add]/t[k])*(1+c*t[k+add]/t[k])))*
//                (-c*(r[k]*t[k+add]/(t[k]*t[k]))/(1+c*t[k+add]/t[k])+c*t[k+add]/(t[k]*t[k])*(r[k+add]+c*r[k]*t[k+add]/t[k])/((1+c*t[k+add]/t[k])*(1+c*t[k+add]/t[k])))*0.25*0.25);
//   }
//
//   double raw_grad_finder(){
//     return (y2[kp]-y1[k])/(x2[kp]-x1[k]);
//   }
//
//   }
//   // double, double, double, double, double, double, double, double, double, double, double, uint, uint, uint
//   void track_coord_true(double, double, double, double, double, double, double, double, double, double, double, uint, uint, uint);
//
//   vector<pair<double,uint> > array_sorter(){
//     vector<pair<double, uint> > vp;
//     for (uint i=0; i<n; i++) {vp.push_back(make_pair(arr[i], i));}
//     sort(vp.begin(), vp.end());
//     return vp;
//   }
//
//   uint * closest_finder(){
//     while (m0<n && m1<n && m2<n && m3<n){
//       double minimum = min(A[m0], min(B[m1], min(C[m2],D[m3])));
//       double maximum = max(A[m0], max(B[m1], max(C[m2],D[m3])));
//
//         if (maximum-minimum < diff){
//           min_idx[0] = m0, min_idx[1] = m1, min_idx[2] = m2, min_idx[3] = m3;
//           diff = maximum-minimum;
//         }
//
//         if (diff == 0) break;
//
//         if      (A[m0] == minimum) m0++;
//         else if (B[m1] == minimum) m1++;
//         else if (C[m2] == minimum) m2++;
//         else                       m3++;
//     }
//     return min_idx;
//
// };
//
// Track::set_track_coord(double r0[], double t0[], uint k0, uint add0, double c0){
//   r[] = r0
// }
//
// Track::track_coord_true(){  if ((n==0) || (n==3)) {
//     x_track[k] = x_track_p[k];
//     y_track[k] = y_track_p[k];
//     x_track_err[k] = track_coord_error(x, t, k, add, 1);
//     y_track_err[k] = track_coord_error(y, t, k, add, 1);
//   }
//   if ((n==1) || (n==2)) {
//     x_track[k] = x_track_m[k];
//     y_track[k] = y_track_m[k];
//     x_track_err[k] = track_coord_error(x, t, k, add, -1);
//     y_track_err[k] = track_coord_error(y, t, k, add, -1);
//   }
}

#endif
