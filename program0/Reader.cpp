#include <bitset>
using namespace std;
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
