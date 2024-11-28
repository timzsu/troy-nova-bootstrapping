#pragma once

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>

#include "RemezArcsin.cuh"
#include "RemezCos.cuh"
#include "common/Polynomial.cuh"

using namespace std;

namespace bootstrap {

class ModularReducer {
 public:
  long boundary_K;
  double log_width;
  long deg;
  long num_double_formula;

  double inverse_log_width;
  long inverse_deg;

  double scale_inverse_coeff;

  CKKSEvaluator *ckks = nullptr;

  RemezParam rmparm;

  RemezCos *poly_generator;
  RemezArcsin *inverse_poly_generator;

  Polynomial sin_cos_polynomial;
  Polynomial inverse_sin_polynomial;

  ModularReducer(long _boundary_K, double _log_width, long _deg, long _num_double_formula, long _inverse_deg, CKKSEvaluator *ckks);

  void double_angle_formula(troy::Ciphertext &cipher);
  void double_angle_formula_scaled(troy::Ciphertext &cipher, double scale_coeff);
  void generate_sin_cos_polynomial();
  void generate_inverse_sine_polynomial();
  void write_polynomials();
  void modular_reduction(troy::Ciphertext &rtn, troy::Ciphertext &cipher);
};

}; // namespace bootstrap