#include <math.h>

double force(double W, double delta, double r){
  return W*delta/(pow(r,3.0));                     //Could seperate into temp=W/pow.. and temp*delta
}



