#pragma once

#include <NTL/RR.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

#include "../ckks_evaluator.cuh"
#include "func.cuh"

using namespace std;
using namespace NTL;

namespace bootstrap {
class Polynomial {
 public:
  RR *coeff = 0;
  long deg = 0, heap_k = 0, heap_m = 0, heaplen = 0;
  RR *chebcoeff = 0;
  Polynomial **poly_heap = 0;

  Polynomial() = default;
  Polynomial(long _deg);
  Polynomial(long _deg, RR *_coeff, string tag);

  ~Polynomial();
  void set_polynomial(long _deg, RR *_coeff, string tag);
  void set_zero_polynomial(long _deg);
  void showcoeff();
  void showchebcoeff();
  void copy(Polynomial &poly);
  void power_to_cheb();
  void cheb_to_power();
  RR evaluate(RR value);

  void constmul(RR constant);
  void mulinplace(Polynomial &poly);
  void addinplace(Polynomial &poly);
  void subtinplace(Polynomial &poly);

  void change_variable_scale(RR scale);  // f_new(x) = f(x / scale)
  void generate_poly_heap_manual(long k, long m);
  void generate_poly_heap();
  void generate_poly_heap_odd();

  void write_heap_to_file(ofstream &out);
  void read_heap_from_file(ifstream &in);

  void homomorphic_poly_evaluation(CKKSEvaluator *ckks, troy::Ciphertext &rtn, troy::Ciphertext &cipher);
};

void mul(Polynomial &rtn, Polynomial &a, Polynomial &b);
void add(Polynomial &rtn, Polynomial &a, Polynomial &b);
void subt(Polynomial &rtn, Polynomial &a, Polynomial &b);

void divide_poly(Polynomial &quotient, Polynomial &remainder, Polynomial &target, Polynomial &divider);
void chebyshev(Polynomial &rtn, long deg);
void second_chebyshev_times_x_for_sine(Polynomial &rtn, long deg);
}  // namespace boot
