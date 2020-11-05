using namespace std;
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
