#pragma once

#include <complex>
#include <memory>
#include <stdexcept>

#include "../ckks_encoder.h"
#include "../evaluator.h"
#include "../encryptor.h"
#include "../decryptor.h"

using namespace std;

namespace bootstrap {

class Encoder {
 private:
  std::shared_ptr<troy::CKKSEncoder> encoder;

 public:
  Encoder() = default;

  Encoder(std::shared_ptr<troy::CKKSEncoder> encoder) {
    this->encoder = encoder;
  }

  inline size_t slot_count() { return encoder->slot_count(); }

  inline void reset_sparse_slots() { 
    throw std::runtime_error("Sparse slot not supported");
    // encoder->reset_sparse_slots(); 
  }
  
  // Vector (of doubles or complexes) inputs
  inline void encode(vector<double> values, troy::ParmsID& parms_id, double scale, troy::Plaintext &plain) {
    if (values.size() == 1) {
      encode(values[0], scale, plain);
      return;
    }
    vector<complex<double>> complex_values(encoder->slot_count(), 0.0 + 0.0i);
    for (size_t i = 0; i < values.size(); i++) {
      complex_values[i] = values[i];
    }
    encoder->encode_complex64_simd(complex_values, parms_id, scale, plain);
  }

  inline void encode(vector<complex<double>> complex_values, double scale, troy::Plaintext &plain) {
    if (complex_values.size() == 1) {
      encode(complex_values[0], scale, plain);
      return;
    }
    complex_values.resize(encoder->slot_count(), 0.0 + 0.0i);
    encoder->encode_complex64_simd(complex_values, std::nullopt, scale, plain);
  }

  inline void encode(vector<double> values, double scale, troy::Plaintext &plain) {
    if (values.size() == 1) {
      encode(values[0], scale, plain);
      return;
    }
    vector<complex<double>> complex_values(encoder->slot_count(), 0.0 + 0.0i);
    for (size_t i = 0; i < values.size(); i++) {
      complex_values[i] = values[i];
    }
    encoder->encode_complex64_simd(complex_values, std::nullopt, scale, plain);
  }


  // Value inputs (fill all slots with that value)
  inline void encode(double value, troy::ParmsID& parms_id, double scale, troy::Plaintext &plain) {
    encoder->encode_float64_single(value, parms_id, scale, plain);
  }

  inline void encode(double value, double scale, troy::Plaintext &plain) {
    encoder->encode_float64_single(value, std::nullopt, scale, plain);
  }

  inline void encode(complex<double> complex_value, double scale, troy::Plaintext &plain) {
    encoder->encode_complex64_single(complex_value, std::nullopt, scale, plain);
  }

  inline void decode(troy::Plaintext &plain, vector<complex<double>> &values) {
    encoder->decode_complex64_simd(plain, values);
  }
  
  inline void decode(troy::Plaintext &plain, vector<double> &values) {
    vector<complex<double>> complex_values(values.size());
    encoder->decode_complex64_simd(plain, complex_values);
    for (size_t i = 0; i < values.size(); i++) {
      values[i] = complex_values[i].real();
    }
  }
};

class Evaluator: public troy::Evaluator {
  private:
    troy::HeContextPointer context;
    Encoder encoder;

  public:

    Evaluator(troy::HeContextPointer context, Encoder encoder) : troy::Evaluator(context) {
      this->context = context;
      this->encoder = encoder;
    }

  // Bootstrapping
  inline void multiply_const(const troy::Ciphertext &ct, double value, troy::Ciphertext &dest) {
    dest = ct;
    multiply_const_inplace(dest, value);
  }

  inline void multiply_const_inplace(troy::Ciphertext &ct, double value) {
    troy::Plaintext const_plain;

    vector<double> values(encoder.slot_count(), value);
    encoder.encode(values, ct.scale(), const_plain);
    mod_switch_plain_to_inplace(const_plain, ct.parms_id());
    multiply_plain_inplace(ct, const_plain);
  }

  inline void add_const(troy::Ciphertext &ct, double value, troy::Ciphertext &dest) {
    dest = ct;
    add_const_inplace(dest, value);
  }

  inline void add_const_inplace(troy::Ciphertext &ct, double value) {
    troy::Plaintext const_plain;

    vector<double> values(encoder.slot_count(), value);
    encoder.encode(values, ct.scale(), const_plain);
    mod_switch_plain_to_inplace(const_plain, ct.parms_id());
    add_plain_inplace(ct, const_plain);
  }

  inline void add_reduced_error(const troy::Ciphertext &ct1, const troy::Ciphertext &ct2, troy::Ciphertext &dest) {
    if (&ct2 == &dest) {
      add_inplace_reduced_error(dest, ct1);
    } else {
      dest = ct1;
      add_inplace_reduced_error(dest, ct2);
    }
  }

  void add_inplace_reduced_error(troy::Ciphertext &ct1, const troy::Ciphertext &ct2);

  inline void sub_reduced_error(const troy::Ciphertext &ct1, const troy::Ciphertext &ct2, troy::Ciphertext &dest) {
    dest = ct1;
    sub_inplace_reduced_error(dest, ct2);
  }

  void sub_inplace_reduced_error(troy::Ciphertext &ct1, const troy::Ciphertext &ct2);

  inline void multiply_reduced_error(const troy::Ciphertext &ct1, const troy::Ciphertext &ct2, const troy::RelinKeys &relin_keys, troy::Ciphertext &dest) {
    if (&ct2 == &dest) {
      multiply_inplace_reduced_error(dest, ct1, relin_keys);
    } else {
      dest = ct1;
      multiply_inplace_reduced_error(dest, ct2, relin_keys);
    }
  }

  void multiply_inplace_reduced_error(troy::Ciphertext &ct1, const troy::Ciphertext &ct2, const troy::RelinKeys &relin_keys);

  inline void double_inplace(troy::Ciphertext &ct) const {
    add_inplace(ct, ct);
  }

  template <typename T, typename = std::enable_if_t<std::is_same<std::remove_cv_t<T>, double>::value || std::is_same<std::remove_cv_t<T>, std::complex<double>>::value>>
  inline void multiply_vector_reduced_error(troy::Ciphertext &ct, std::vector<T> &values, troy::Ciphertext &dest) {
    dest = ct;
    multiply_vector_inplace_reduced_error(dest, values);
  }

  inline void multiply_vector_inplace_reduced_error(troy::Ciphertext &ct, vector<double> &values) {
    troy::Plaintext plain;

    values.resize(encoder.slot_count(), 0.0);
    encoder.encode(values, ct.scale(), plain);
    mod_switch_plain_to_inplace(plain, ct.parms_id());
    multiply_plain_inplace(ct, plain);
  }

  inline void multiply_vector_inplace_reduced_error(troy::Ciphertext &ct, vector<complex<double>> &values) {
    troy::Plaintext plain;

    values.resize(encoder.slot_count(), 0.0 + 0.0i);
    encoder.encode(values, ct.scale(), plain);
    mod_switch_plain_to_inplace(plain, ct.parms_id());
    multiply_plain_inplace(ct, plain);
  }
};

class CKKSEvaluator {
 private:
  // Sign function coefficients
  vector<double> F4_COEFFS = {0, 315, 0, -420, 0, 378, 0, -180, 0, 35};
  int F4_SCALE = (1 << 7);
  vector<double> G4_COEFFS = {0, 5850, 0, -34974, 0, 97015, 0, -113492, 0, 46623};
  int G4_SCALE = (1 << 10);

  // Helper functions
  uint64_t get_modulus(troy::Ciphertext &x, int k);

  troy::Ciphertext init_guess(troy::Ciphertext x);
  troy::Ciphertext eval_line(troy::Ciphertext x, troy::Plaintext m, troy::Plaintext c);

  // Evaluation functions
  troy::Ciphertext newton_iter(troy::Ciphertext x, troy::Ciphertext res, int iter);
  pair<troy::Ciphertext, troy::Ciphertext> goldschmidt_iter(troy::Ciphertext v, troy::Ciphertext y, int d = 1);
  void eval_odd_deg9_poly(vector<double> &a, troy::Ciphertext &x, troy::Ciphertext &dest);

 public:
  // Memory managed outside of the evaluator
  troy::RelinKeys relin_keys;
  troy::GaloisKeys galois_keys;
  std::vector<std::uint32_t> galois_elts;

  // Component classes
  troy::HeContextPointer context;
  Encoder encoder;
  std::shared_ptr<troy::Encryptor> encryptor;
  std::shared_ptr<Evaluator> evaluator;
  std::shared_ptr<troy::Decryptor> decryptor;

  size_t degree;
  double scale;
  size_t slot_count;

  CKKSEvaluator(troy::EncryptionParameters *parms, double scale, vector<uint32_t> galois_elts = {}) {
    context = troy::HeContext::create(*parms, true, troy::SecurityLevel::Classical128);

    this->scale = scale;
    troy::CKKSEncoder tmp(context); 
    auto encoder_internal = std::make_shared<troy::CKKSEncoder>(std::move(tmp));
    this->slot_count = encoder_internal->slot_count();
    this->degree = this->slot_count * 2;
    if (troy::utils::device_count() > 0) {
        context->to_device_inplace();
        encoder_internal->to_device_inplace();
    }

    troy::KeyGenerator keygen(context);
    auto secret_key = keygen.secret_key();
    troy::PublicKey public_key = keygen.create_public_key(false);
    this->relin_keys = keygen.create_relin_keys(false);
    this->galois_keys = keygen.create_galois_keys(false);
    this->galois_elts = galois_elts;

    // Instantiate the component classes
    this->encoder = Encoder(encoder_internal);
    this->encryptor = std::make_shared<troy::Encryptor>(troy::Encryptor(context)); 
    this->encryptor->set_public_key(public_key);
    this->evaluator = std::make_shared<Evaluator>(Evaluator(context, encoder)); 
    this->decryptor = std::make_shared<troy::Decryptor>(troy::Decryptor(context, secret_key));
  }

  // Helper functions
  vector<double> init_vec_with_value(double value);
  troy::Plaintext init_plain_power_of_x(size_t exponent);

  void re_encrypt(troy::Ciphertext &ct);
  void print_decrypted_ct(troy::Ciphertext &ct, int num);
  void print_decoded_pt(troy::Plaintext &pt, int num);

  // Evaluation functions
  troy::Ciphertext sgn_eval(troy::Ciphertext x, int d_g, int d_f, double sgn_factor = 0.5);
  troy::Ciphertext invert_sqrt(troy::Ciphertext x, int d_newt = 20, int d_gold = 1);
  troy::Ciphertext exp(troy::Ciphertext x);
  troy::Ciphertext inverse(troy::Ciphertext x, int iter = 4);

  // Metrics calcuation functions
  double calculate_MAE(vector<double> &y_true, troy::Ciphertext &ct, int N);
};

}; // namespace bootstrap