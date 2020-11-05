#include <iostream>
#include <bitset>
#include <cmath>
#include <chrono>
#include <tuple>
#include <vector>
#include <cfloat>
#include <algorithm>

#include "Reader.cpp"
#include "TrackCoord.cpp"
#include "Gradient.cpp"
#include "Fitter.cpp"
#include "Parameters.cpp"
#include "OneEvent.cpp"
// #include <Ploter.cpp>
using namespace std;
using namespace std::chrono;
typedef unsigned int uint;


// ******************************************************************************
// main
// ******************************************************************************
int main(int argc, char *argv[]) {

  // open the data file
  FILE * f=fopen(("%s",argv[1]),"rb");

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
    // hist(vel_new,   "c1", "hist1", "Velocity [ #mum/ns ]", "drift_velocity.pdf",    "#mum/ns",  51.5, 53.5);
    // hist(angle_new, "c2", "hist2", "Angle [ #circ ]",      "inclination_angle.pdf", "#circ",   -15,   25);
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
