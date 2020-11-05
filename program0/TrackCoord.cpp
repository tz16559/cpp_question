

using namespace std;
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
