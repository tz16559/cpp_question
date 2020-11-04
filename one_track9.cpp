#include <iostream>
#include <bitset>
#include <cmath>
#include <chrono>
#include <tuple>
#include <bits/stdc++.h>
using namespace std;
using namespace std::chrono;
typedef unsigned int uint;


double track_coord_finder(double r[], double t[], uint k, uint add, double c){
  return (r[k+add]+c*r[k]*t[k+add]/t[k])/(1+c*t[k+add]/t[k]);
}

double track_coord_error(double r[], double t[], uint k, uint add, double c){
  return sqrt((c*(r[k]/t[k])/(1+c*t[k+add]/t[k])-c*1/t[k]*(r[k+add]+c*r[k]*t[k+add]/t[k])/((1+c*t[k+add]/t[k])*(1+c*t[k+add]/t[k])))*
              (c*(r[k]/t[k])/(1+c*t[k+add]/t[k])-c*1/t[k]*(r[k+add]+c*r[k]*t[k+add]/t[k])/((1+c*t[k+add]/t[k])*(1+c*t[k+add]/t[k])))*0.25*0.25
            +(-c*(r[k]*t[k+add]/(t[k]*t[k]))/(1+c*t[k+add]/t[k])+c*t[k+add]/(t[k]*t[k])*(r[k+add]+c*r[k]*t[k+add]/t[k])/((1+c*t[k+add]/t[k])*(1+c*t[k+add]/t[k])))*
             (-c*(r[k]*t[k+add]/(t[k]*t[k]))/(1+c*t[k+add]/t[k])+c*t[k+add]/(t[k]*t[k])*(r[k+add]+c*r[k]*t[k+add]/t[k])/((1+c*t[k+add]/t[k])*(1+c*t[k+add]/t[k])))*0.25*0.25);
}

double raw_grad_finder(double x1[], double y1[], double x2[], double y2[], uint k, uint kp){
  return (y2[kp]-y1[k])/(x2[kp]-x1[k]);
}

tuple<double, double, double, double> line_fitter(double x[], double y[], double x_err[], double y_err[], uint n){
  double coef,cons,coef_err,cons_err;
  double s_x = 0, s_y=0, s_xy=0, s_xx=0, s=0,s_tt=0;
  double *a, *sigma;
  a = new double [n]; sigma = new double [n];

  for (uint i=0; i<n; i++){sigma[i] = sqrt(x_err[i]*x_err[i]+y_err[i]*y_err[i]);}
  for (uint i=0; i<n; i++){
    s_xy += x[i]*y[i]/(sigma[i]*sigma[i]);
    s_xx += x[i]*x[i]/(sigma[i]*sigma[i]);
    s_x  += x[i]     /(sigma[i]*sigma[i]);
    s_y  += y[i]     /(sigma[i]*sigma[i]);
    s    += 1        /(sigma[i]*sigma[i]);
  }
  for (uint i=0; i<n; i++){a[i] = 1/sigma[i]*(x[i]-s_x/s);}
  for (uint i=0; i<n; i++){s_tt += a[i]*a[i];}

  for (uint i=0; i<n; i++){coef += 1/s_tt*(a[i]*y[i]/sigma[i]);}//m
                           cons = (s_y-s_x*coef)/s; //c
                           coef_err = sqrt(1/s_tt);
                           cons_err = sqrt(1/s*(1+s_x*s_x/(s*s_tt)));
  delete[] a,sigma;
  return make_tuple(coef,cons,coef_err,cons_err);
}

tuple<int,int,int,int> closest_finder(double A[], double B[], double C[], double D[], uint n){
  double diff=DBL_MAX;
  uint m0=0, m1=0, m2=0, m3=0;
  uint min0=0, min1=0, min2=0, min3=0;
  // printf("eo");
  while (m0<n && m1<n && m2<n && m3<n){
    double minimum = min(A[m0], min(B[m1], min(C[m2],D[m3])));
    double maximum = max(A[m0], max(B[m1], max(C[m2],D[m3])));
    // printf("%f",minimum);
      if (maximum-minimum < diff){
        min0 = m0, min1 = m1, min2 = m2, min3 = m3;
        diff = maximum-minimum;
      }

       if (diff == 0) break;
      //{
      //   min0 =0, min1 = 0, min2 = 0, min3 =0;
      // }

      if      (A[m0] == minimum) m0++;
      else if (B[m1] == minimum) m1++;
      else if (C[m2] == minimum) m2++;
      else                       m3++;
  }

  return make_tuple(min0, min1, min2, min3);
}

void track_coord_true(double x_track[],   double y_track[],   double x_track_err[], double y_track_err[],
                      double x_track_p[], double y_track_p[], double x_track_m[],   double y_track_m[],
                      double t[],         double x[],         double y[],
                      uint n,             uint k,             uint add){

  if ((n==0) || (n==3)) {
    x_track[k] = x_track_p[k];
    y_track[k] = y_track_p[k];
    x_track_err[k] = track_coord_error(x, t, k, add, 1);
    y_track_err[k] = track_coord_error(y, t, k, add, 1);
  }
  if ((n==1) || (n==2)) {
    x_track[k] = x_track_m[k];
    y_track[k] = y_track_m[k];
    x_track_err[k] = track_coord_error(x, t, k, add, -1);
    y_track_err[k] = track_coord_error(y, t, k, add, -1);
  }
}

vector<pair<double,uint> > array_sorter(double arr[], uint n){
  vector<pair<double, uint> > vp;
  for (uint i=0; i<n; i++) {vp.push_back(make_pair(arr[i], i));}
  sort(vp.begin(), vp.end());
  vp.clear();
  return vp;
}

tuple<double, double> incl_angle_finder(double m, double m_err){
  double angle     = 180/M_PI*atan(m);
  double angle_err = 180/M_PI/(1+m*m)*m_err;
  return make_tuple(angle, angle_err);
}

tuple<double, double> avg_vel_finder(double m, double m_err, double c, double c_err, double t[], double x[], double y[]){
  double *vel, *vel_err, *vel_wgt;
  vel = new double [8]; vel_err= new double [8]; vel_wgt= new double [8];
  double sum_vel_wgt=0, vel_avg=0, vel_avg_err;

  for (uint i=0; i<8; i++){
    vel[i]     = 1/(t[i])*abs(m*x[i]-y[i]+c)/(sqrt(m*m+1));
    vel_err[i] = sqrt((vel[i]*0.25/t[i])*
                      (vel[i]*0.25/t[i])
                       +(((m*x[i]-y[i]+c)*x[i])/(vel[i]*t[i]*t[i]*(m*m+1))-(vel[i]*m)/(m*m+1))*
                        (((m*x[i]-y[i]+c)*x[i])/(vel[i]*t[i]*t[i]*(m*m+1))-(vel[i]*m)/(m*m+1))*m_err*m_err
                           +((m*x[i]-y[i]+c)/(t[i]*t[i]*vel[i]*(m*m+1)))*
                            ((m*x[i]-y[i]+c)/(t[i]*t[i]*vel[i]*(m*m+1)))*c_err*c_err);
    vel_wgt[i] = 1/vel_err[i];
  }

  for (uint i=0; i<8; i++){sum_vel_wgt += vel_wgt[i];}
  for (uint i=0; i<8; i++){vel_avg  += (vel_wgt[i]*vel[i])/sum_vel_wgt;}
  vel_avg_err = 1/(sqrt(sum_vel_wgt));
  delete[] vel, vel_err, vel_wgt;
  return make_tuple(vel_avg, vel_avg_err);
}

tuple<double, double, double, double> run_one_event(FILE * f, uint event){
  uint16_t *hits; hits = new uint16_t[8];
  fseek(f, event*sizeof(uint16_t)*8, SEEK_SET);
  fread(hits, 2, 8, f);
  double *t, *x, *y; t = new double[8]; x= new double[8]; y = new double[8];
  for (uint i=0; i<8; i++){
    bitset<16> binary(hits[i]);
    t[i] = (hits[i] >> 6) * 0.5 + 0.25;
    x[i] = 0b111 & hits[i];
    y[i] = (0b111000 & hits[i]) >> 3;
    int x_int = int(x[i]);
    if (x_int % 2 != 0) y[i] += 0.5;
    x[i] *= pow(10,4);
    y[i] *= pow(10,4);
    // cout << t[i] << x[i]<< y[i] <<event <<endl;
  }

  uint   n=4, m=4, add=4;
  double *x_track_p, *y_track_p, *x_track_m, *y_track_m;
  x_track_p = new double[n];
  y_track_p = new double[n];
  x_track_m = new double[n];
  y_track_m = new double[n];
  double **grad = new double *[n];
  for (uint i=0; i<n;i++){
    grad[i] = new double [m];
  }


  // double x_track_p[n], y_track_p[n];
  // double x_track_m[n], y_track_m[n];
  // double grad[n][m];
  //
  for (uint k=0; k<n; k++){
    x_track_p[k] = track_coord_finder(x, t, k, add, 1);
    y_track_p[k] = track_coord_finder(y, t, k, add, 1);
    if (t[k] != t[k+add]){
      x_track_m[k] = track_coord_finder(x, t, k, add, -1);
      y_track_m[k] = track_coord_finder(y, t, k, add, -1);
    }
    else{
      x_track_m[k] = x_track_p[k];
      y_track_m[k] = y_track_p[k];
    }
  }

  for (uint k=0; k<n; k++){
    uint kp = k+1;
    if ((kp) >= n) kp -= 4;
    grad[k][0] = raw_grad_finder(x_track_p, y_track_p, x_track_p, y_track_p, k, kp);
    grad[k][1] = raw_grad_finder(x_track_m, y_track_m, x_track_m, y_track_m, k, kp);
    grad[k][2] = raw_grad_finder(x_track_m, y_track_m, x_track_p, y_track_p, k, kp);
    grad[k][3] = raw_grad_finder(x_track_p, y_track_p, x_track_m, y_track_m, k, kp);
    // printf("%f %f %f %f\n", grad[k][0], grad[k][1], grad[k][2], grad[k][3]);
  }


  vector<pair<double, uint> > g0_sort = array_sorter(grad[0], m);
  vector<pair<double, uint> > g1_sort = array_sorter(grad[1], m);
  vector<pair<double, uint> > g2_sort = array_sorter(grad[2], m);
  vector<pair<double, uint> > g3_sort = array_sorter(grad[3], m);


  double *grad0_sorted, *grad0_sorted_index;
  double *grad1_sorted, *grad1_sorted_index;
  double *grad2_sorted, *grad2_sorted_index;
  double *grad3_sorted, *grad3_sorted_index;

  grad0_sorted = new double [m]; grad0_sorted_index = new double [m];
  grad1_sorted = new double [m]; grad1_sorted_index = new double [m];
  grad2_sorted = new double [m]; grad2_sorted_index = new double [m];
  grad3_sorted = new double [m]; grad3_sorted_index = new double [m];

  for (uint g=0; g<m; g++){
    grad0_sorted[g] = g0_sort[g].first; grad0_sorted_index[g] = g0_sort[g].second;
    grad1_sorted[g] = g1_sort[g].first; grad1_sorted_index[g] = g1_sort[g].second;
    grad2_sorted[g] = g2_sort[g].first; grad2_sorted_index[g] = g2_sort[g].second;
    grad3_sorted[g] = g3_sort[g].first; grad3_sorted_index[g] = g3_sort[g].second;
  }

  uint gs0, gs1, gs2, gs3, g0, g1, g2, g3;
  tie(gs0, gs1, gs2, gs3) = closest_finder(grad0_sorted, grad1_sorted, grad2_sorted, grad3_sorted, m);
  // return make_tuple (grad[0][0], grad[0][1], grad[0][2], grad[0][3]);

  for (uint g=0; g<4; g++){
    if (gs0 == g) {g0 = grad0_sorted_index[g];}
    if (gs0 == g) {g1 = grad1_sorted_index[g];}
    if (gs0 == g) {g2 = grad2_sorted_index[g];}
    if (gs0 == g) {g3 = grad3_sorted_index[g];}
  }

  double *grad_true, *x_track, *y_track, *x_track_err, *y_track_err;
  grad_true = new double[n];
  x_track = new double[n];
  y_track = new double[n];
  x_track_err = new double[n];
  y_track_err = new double[n];

  grad_true[0] = grad[0][g0];
  grad_true[1] = grad[1][g1];
  grad_true[2] = grad[2][g2];
  grad_true[3] = grad[3][g3];

  track_coord_true(x_track, y_track, x_track_err, y_track_err, x_track_p, y_track_p, x_track_m, y_track_m, t, x, y, g0, 0, add);
  track_coord_true(x_track, y_track, x_track_err, y_track_err, x_track_p, y_track_p, x_track_m, y_track_m, t, x, y, g1, 1, add);
  track_coord_true(x_track, y_track, x_track_err, y_track_err, x_track_p, y_track_p, x_track_m, y_track_m, t, x, y, g2, 2, add);
  track_coord_true(x_track, y_track, x_track_err, y_track_err, x_track_p, y_track_p, x_track_m, y_track_m, t, x, y, g3, 3, add);

  double coef,cons,coef_err,cons_err, angle, angle_err, vel_avg, vel_avg_err;
  tie(coef,cons,coef_err,cons_err) = line_fitter(x_track, y_track, x_track_err, y_track_err, n);
  tie(angle,angle_err)  = incl_angle_finder(coef, coef_err);
  tie(vel_avg,vel_avg_err) = avg_vel_finder(coef, coef_err, cons, cons_err, t, x, y);

  for(uint i=0; i<n; ++i) {
    delete[] grad[i];
}
// Free the array of pointers
delete[] hits, t, x, y;
delete[] grad, x_track_p, y_track_p, x_track_m, y_track_m,
         grad0_sorted, grad0_sorted_index, grad1_sorted, grad1_sorted_index, grad2_sorted, grad2_sorted_index, grad3_sorted, grad3_sorted_index,
         grad_true, x_track, y_track, x_track_err, y_track_err;

  // printf("%f %f %f %f", coef, coef_err, cons, cons_err);
  return make_tuple(angle, angle_err, vel_avg, vel_avg_err);
}

int main() {
  // FILE * f=fopen("onetrack.raw","rb");

  uint e = 1000000;
  double *angle, *angle_err, *vel, *vel_err;

  // auto start = high_resolution_clock::now();
  double angle0, angle_err0, vel0, vel_err0,angle1, angle_err1, vel1, vel_err1, angle2, angle_err2, vel2, vel_err2, angle3, angle_err3, vel3, vel_err3;
  FILE * f=fopen("manytracks.raw","rb");
  tie(angle0, angle_err0, vel0, vel_err0) = run_one_event(f,999999);
  angle = new double [e], angle_err= new double [e], vel= new double [e], vel_err= new double [e];
  printf("%d angle = %f(%f), mean_velelocity = %f(%f)\n", 0, angle0, angle_err0, vel0, vel_err0);
  uint i=0;
  while (i<e){


    tie(angle[i], angle_err[i], vel[i], vel_err[i]) = run_one_event(f,i);
    // output[i]  = run_one_event(f,i);

     printf("%d angle = %f(%f), mean_velelocity = %f(%f)\n", i, angle[i], angle_err[i], vel[i], vel_err[i]);
     i +=1;
  }

  fclose(f);
  // auto stop = high_resolution_clock::now();
  // auto duration = duration_cast<microseconds>(stop - start);
  // cout << "Time taken by function: "
  //      << duration.count() << " microseconds" << endl;

  return 0;
}



// //g++ -o root.exe root_test.cpp $(root-config --cflags --libs)
//
