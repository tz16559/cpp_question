using namespace std;
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
