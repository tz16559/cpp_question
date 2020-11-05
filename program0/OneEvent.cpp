using namespace std;
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
