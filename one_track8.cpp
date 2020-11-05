#include <iostream>
#include <bitset>
#include <cmath>
#include <chrono>
#include <tuple>
#include <vector>

//Plotting
#include "TH1.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TGaxis.h"

using namespace std;
using namespace std::chrono;
typedef unsigned int uint;

void event_reader(FILE *f, uint event, uint h, double t[], double x[], double y[], uint16_t hits[]){
  // read one event containing h=8 hits from the whole file
  fseek(f, event*sizeof(uint16_t)*h, SEEK_SET);
  fread(hits, sizeof(uint16_t), h, f);
  for (uint i=0; i<h; i++){
    bitset<16> binary(hits[i]);
    t[i] = (hits[i] >> 6) * 0.5 + 0.25;
    x[i] = 0b111 & hits[i];
    y[i] = (0b111000 & hits[i]) >> 3;
    int x_int = int(x[i]);
    if (x_int % 2 != 0) y[i] += 0.5;
    x[i] *= pow(10,4);
    y[i] *= pow(10,4);
  }
}

double track_coord_finder(double r[], double t[], uint k, uint add, double c){
  // assign values to the x or y coordinate of one point on the track
  // use the tangential ratios between the addresses of two wires
  return (r[k+add]+c*r[k]*t[k+add]/t[k])/(1+c*t[k+add]/t[k]);
}

void track_coord_all(double x_track_p[], double y_track_p[], double x_track_m[], double y_track_m[],
                     double x[],         double y[],         double t[],
                     uint n,             uint add){
  // assign values of coordinates for all possible points on the track: (xp,yp) and (xm,ym)
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
}

double track_coord_error(double r[], double t[], uint k, uint add, double c){
  // calculate the uncertainties in the coordinates of the points on the track
  return sqrt((c*(r[k]/t[k])/(1+c*t[k+add]/t[k])-c*1/t[k]*(r[k+add]+c*r[k]*t[k+add]/t[k])/((1+c*t[k+add]/t[k])*(1+c*t[k+add]/t[k])))*
              (c*(r[k]/t[k])/(1+c*t[k+add]/t[k])-c*1/t[k]*(r[k+add]+c*r[k]*t[k+add]/t[k])/((1+c*t[k+add]/t[k])*(1+c*t[k+add]/t[k])))*0.25*0.25
            +(-c*(r[k]*t[k+add]/(t[k]*t[k]))/(1+c*t[k+add]/t[k])+c*t[k+add]/(t[k]*t[k])*(r[k+add]+c*r[k]*t[k+add]/t[k])/((1+c*t[k+add]/t[k])*(1+c*t[k+add]/t[k])))*
             (-c*(r[k]*t[k+add]/(t[k]*t[k]))/(1+c*t[k+add]/t[k])+c*t[k+add]/(t[k]*t[k])*(r[k+add]+c*r[k]*t[k+add]/t[k])/((1+c*t[k+add]/t[k])*(1+c*t[k+add]/t[k])))*0.25*0.25);
}

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

void track_coord_true(double x_track[],   double y_track[],   double x_track_err[], double y_track_err[],
                      double x_track_p[], double y_track_p[], double x_track_m[],   double y_track_m[],
                      double t[],         double x[],         double y[],
                      uint n,             uint k,             uint add){
// find the coordinates on the track that gives the combinations to have the closest values in gradient
// this is the coordinates we are intersted in
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

tuple<double, double, double, double> line_fitter(double x[], double y[], double x_err[], double y_err[], uint n){
  // least squares fit for a straight line y=mx+c
  // find the equation of the track using the true coordinates found above
  double coef,cons,coef_err,cons_err;
  double s_x = 0, s_y=0, s_xy=0, s_xx=0, s=0, s_tt=0;
  double a[n], sigma[n];
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
  return make_tuple(coef,cons,coef_err,cons_err);
}

double * incl_angle_finder(double m, double m_err){
  // find the inclination angle of the track from the horizontal
  static double angle[2];
  angle[0] = 180/M_PI*atan(m);
  angle[1] = 180/M_PI/(1+m*m)*m_err;
  return angle;
}

double * avg_vel_finder(double m, double m_err, double c, double c_err, double t[], double x[], double y[]){
  // find the averaged drift velocity for each track
  // each velocity was calculated from the tangential distance between the track and the each wire with the corresponding drift time
  double vel[8], vel_err[8], vel_wgt[8], sum_vel_wgt=0;
  static double vel_avg[2]; vel_avg[0]=0;

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
  for (uint i=0; i<8; i++){vel_avg[0]  += (vel_wgt[i]*vel[i])/sum_vel_wgt;}
  vel_avg[1] = 1/(sqrt(sum_vel_wgt));
  return vel_avg;
}

void hist(vector<double> var,
          const char* c_no,  const char *hist_no,
          const char *title, const char *name, const char *units,
          double low_b,      double high_b){
  // plot the results into histograms
  TGaxis::SetMaxDigits(4);
  TCanvas *c1 = new TCanvas(c_no, "c1", 150, 10, 990, 660);
  TH1F *hist  = new TH1F(hist_no, " ", 200, low_b, high_b);
  for (uint i=0;i<var.size();i++) {hist->Fill(var[i]);}

  hist->GetXaxis()->SetTitle(title);
  hist->GetYaxis()->SetTitle(Form("Candidates/(%.2f %s)",hist->GetBinWidth(1),units));
  hist->Draw("E");
  c1->SaveAs(name);
}

tuple<double, double, double, double> one_event_analyser(FILE * f, uint event){
// analysis for one event

// ****************************************************************************
// declare the variables
// ****************************************************************************

  // index variables
  uint     h=8;
  uint     n=4;
  uint     m=4;
  uint     add=4;

  // hits
  uint16_t hits[h];

  // drift_time and wire addresses
  double t[h];
  double x[h];
  double y[h];

  // track coordinates
  double x_track_p[n];
  double y_track_p[n];
  double x_track_m[n];
  double y_track_m[n];
  double x_track[n];
  double y_track[n];
  double x_track_err[n];
  double y_track_err[n];

  // all possible gradients, and variables/indices used to find the correct gradients
  double grad[n][m];
  double grad_true[n];

  uint *gs;
  uint g0, g1, g2, g3;
  vector<pair<double,uint> > g0_sort;
  vector<pair<double,uint> > g1_sort;
  vector<pair<double,uint> > g2_sort;
  vector<pair<double,uint> > g3_sort;

  double grad0_sorted[m];
  double grad1_sorted[m];
  double grad2_sorted[m];
  double grad3_sorted[m];
  double grad0_sorted_index[m];
  double grad1_sorted_index[m];
  double grad2_sorted_index[m];
  double grad3_sorted_index[m];

  // final results in fitting the equation of the track, drift velocity and the inclination angle
  double coef, cons, coef_err, cons_err;
  double *angle;
  double *vel_avg;

// ****************************************************************************
// start the analysis
// ****************************************************************************

  // read the events and find the possible coordinates on the track
  event_reader(f, event, h, t, x, y, hits);
  track_coord_all(x_track_p, y_track_p, x_track_m, y_track_m, x, y, t, n, add);

  // find the possible gradients for the coordinates combinations
  for (uint k=0; k<n; k++){
    uint kp = k+1;
    if ((kp) >= n) kp -= 4;
    grad[k][0] = raw_grad_finder(x_track_p, y_track_p, x_track_p, y_track_p, k, kp);
    grad[k][1] = raw_grad_finder(x_track_m, y_track_m, x_track_m, y_track_m, k, kp);
    grad[k][2] = raw_grad_finder(x_track_m, y_track_m, x_track_p, y_track_p, k, kp);
    grad[k][3] = raw_grad_finder(x_track_p, y_track_p, x_track_m, y_track_m, k, kp);
  }

  // find the closest value from each set of gradients by sorting the set first
  // this gives the corrent gradient and thus the correct combination of coordinates
  g0_sort = array_sorter(grad[0], m);
  g1_sort = array_sorter(grad[1], m);
  g2_sort = array_sorter(grad[2], m);
  g3_sort = array_sorter(grad[3], m);

  for (uint g=0; g<m; g++){
    grad0_sorted[g] = g0_sort[g].first; grad0_sorted_index[g] = g0_sort[g].second;
    grad1_sorted[g] = g1_sort[g].first; grad1_sorted_index[g] = g1_sort[g].second;
    grad2_sorted[g] = g2_sort[g].first; grad2_sorted_index[g] = g2_sort[g].second;
    grad3_sorted[g] = g3_sort[g].first; grad3_sorted_index[g] = g3_sort[g].second;
  }

  gs = closest_finder(grad0_sorted, grad1_sorted, grad2_sorted, grad3_sorted, m);

  for (uint g=0; g<4; g++){
    if (gs[0] == g) {g0 = grad0_sorted_index[g];}
    if (gs[1] == g) {g1 = grad1_sorted_index[g];}
    if (gs[2] == g) {g2 = grad2_sorted_index[g];}
    if (gs[3] == g) {g3 = grad3_sorted_index[g];}
  }

  grad_true[0] = grad[0][g0];
  grad_true[1] = grad[1][g1];
  grad_true[2] = grad[2][g2];
  grad_true[3] = grad[3][g3];

  track_coord_true(x_track, y_track, x_track_err, y_track_err, x_track_p, y_track_p, x_track_m, y_track_m, t, x, y, g0, 0, add);
  track_coord_true(x_track, y_track, x_track_err, y_track_err, x_track_p, y_track_p, x_track_m, y_track_m, t, x, y, g1, 1, add);
  track_coord_true(x_track, y_track, x_track_err, y_track_err, x_track_p, y_track_p, x_track_m, y_track_m, t, x, y, g2, 2, add);
  track_coord_true(x_track, y_track, x_track_err, y_track_err, x_track_p, y_track_p, x_track_m, y_track_m, t, x, y, g3, 3, add);

  // calculate the final results
  tie(coef, cons, coef_err, cons_err) = line_fitter(x_track, y_track, x_track_err, y_track_err, n);
  angle   = incl_angle_finder(coef, coef_err);
  vel_avg = avg_vel_finder(coef, coef_err, cons, cons_err, t, x, y);

  return make_tuple(angle[0], angle[1], vel_avg[0], vel_avg[1]);
}

// ******************************************************************************
// main
// ******************************************************************************
int main(int argc, char *argv[]) {

  // open the data file
  FILE * f=fopen(Form("%s",argv[1]),"rb");

  // store the results into arrays
  double angle0;
  double angle_err0;
  double vel0;
  double vel_err0;

  double *angle;
  double *angle_err;
  double *vel;
  double *vel_err;
  uint e    = 1000000;
  angle     = new double [e];
  angle_err = new double [e];
  vel       = new double [e];
  vel_err   = new double [e];

  // results after rejecting problematic events
  vector<double> angle_new;
  vector<double> angle_new_err;
  vector<double> vel_new;
  vector<double> vel_new_err;

  // start the timer
  auto start = high_resolution_clock::now();

  // make the analysis for different data files
  if (string(argv[1]) == "onetrack.raw"){
    tie(angle0, angle_err0, vel0, vel_err0)=one_event_analyser(f, 0);
    printf("Results for one track:\nmean drift velocity = %f +- %f um/ns,\ninclination angle = %f +- %f degree\n", angle0, angle_err0, vel0, vel_err0);
  }

  if (string(argv[1]) == "manytracks.raw"){
    printf("Results for %d tracks, analysing ...\n", e);
    for (uint i=0; i<e; i++){
      tie(angle[i], angle_err[i], vel[i], vel_err[i]) = one_event_analyser(f,i);
      uint j=0;
      if (angle[i]!=angle[i]) j=1;
      if (j!=1){
        angle_new.push_back(angle[i]);
        angle_new_err.push_back(angle_err[i]);
        vel_new.push_back(vel[i]);
        vel_new_err.push_back(vel_err[i]);
      }
    }

    // generate plots for the distributons of the drift velocity and the angle
    hist(vel_new,   "c1", "hist1", "Velocity [ #mum/ns ]", "drift_velocity.pdf",    "#mum/ns",  51.5, 53.5);
    hist(angle_new, "c2", "hist2", "Angle [ #circ ]",      "inclination_angle.pdf", "#circ",   -15,   25);
  }

  // stop the timer
  auto stop = high_resolution_clock::now();

  // close the data file
  fclose(f);

  // print out the execute time for the whole program
  auto duration = duration_cast<microseconds>(stop - start);
  cout << "Execute time: "
       << duration.count()*pow(10,-6) << " seconds" << endl;

  return 0;
}

//
//  g++ -o output8 one_track8.cpp $(root-config --cflags --libs)
//  ./output8 onetrack.raw
//  ./output8 manytracks.raw
