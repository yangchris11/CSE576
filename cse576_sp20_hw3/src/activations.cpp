
#include <cmath>
#include <algorithm>

#include "matrix.h"
#include "neural.h"

#include <iostream>

// Calculate a linear activation (i.e. no activation).
//  f(x) = x
// Parameters:
//  Matrix &matrix: the input non-activated output of the layer.
// Returns:
//  A Matrix containing the activated output.
Matrix forward_linear(const Matrix &matrix) {
  Matrix activated = matrix;
  // TODO: Implement forward activation.
  return activated;
}

// Calculate the backwards pass for the activation.
// Parameters:
//  const Matrix& out: the activated output of the current layer.
//  const Matrix& prev_grad: the gradient from the next layer (towards the Loss).
// Returns:
//  Matrix: the gradients of this layer (to be passed to the previous layer).
Matrix backward_linear(const Matrix &out, const Matrix &prev_grad) {
  assert_same_size(prev_grad, out);
  Matrix grad = prev_grad;
  // TODO: Implement activation backward pass.
  return grad;
}

// Calculate a logistic activation (sigmoid).
// Parameters:
//  Matrix &matrix: the input non-activated output of the layer.
// Returns:
//  A Matrix containing the activated output.
Matrix forward_logistic(const Matrix &matrix) {
  Matrix activated = matrix;
  // TODO: Implement forward activation.
  // Hint: look at matrix.h, it might save you some typing.
  for(int i = 0 ; i < matrix.rows ; ++i){
    for(int j = 0 ; j < matrix.cols ; ++j){
      activated(i, j) = 1.0 / (1 + std::exp(-1 * matrix(i, j)));
    }
  }
  return activated;
}

// Calculate the backwards pass for the activation.
// Parameters:
//  const Matrix& out: the activated output of the current layer.
//  const Matrix& prev_grad: the gradient from the next layer (towards the Loss).
// Returns:
//  Matrix: the gradients of this layer (to be passed to the previous layer).
Matrix backward_logistic(const Matrix &out, const Matrix &prev_grad) {
  assert_same_size(prev_grad, out);
  Matrix grad = prev_grad;
  // TODO: Implement activation backward pass.
  for(int i = 0 ; i < grad.rows ; ++i){
    for(int j = 0 ; j < grad.cols ; ++j){
      grad(i, j) = prev_grad(i, j) * (out(i, j) * (1 - out(i, j)));
    }
  }
  return grad;
}

// Calculate a tanh activation.
// Parameters:
//  Matrix &matrix: the input non-activated output of the layer.
// Returns:
//  A Matrix containing the activated output.
Matrix forward_tanh(const Matrix &matrix) {
  Matrix activated = matrix;
  // TODO: Implement forward activation.
  for(int i = 0 ; i < matrix.rows ; ++i){
    for(int j = 0 ; j < matrix.cols ; ++j){
      activated(i, j) = std::tanh(matrix(i, j));
    }
  }
  return activated;
}

// Calculate the backwards pass for the activation.
// Parameters:
//  const Matrix& out: the activated output of the current layer.
//  const Matrix& prev_grad: the gradient from the next layer (towards the Loss).
// Returns:
//  Matrix: the gradients of this layer (to be passed to the previous layer).
Matrix backward_tanh(const Matrix &out, const Matrix &prev_grad) {
  assert_same_size(prev_grad, out);
  Matrix grad = prev_grad;
  // TODO: Implement activation backward pass.
  for(int i = 0 ; i < grad.rows ; ++i){
    for(int j = 0 ; j < grad.cols ; ++j){
      grad(i, j) = prev_grad(i, j) * (1 - std::pow(out(i, j), 2.0));
    }
  }
  return grad;
}

// Calculate a ReLU activation.
// Parameters:
//  Matrix &matrix: the input non-activated output of the layer.
// Returns:
//  A Matrix containing the activated output.
Matrix forward_relu(const Matrix &matrix) {
  // TODO: Implement forward activation.
  Matrix activated = matrix;
  for(int i = 0 ; i < matrix.rows ; ++i){
    for(int j = 0 ; j < matrix.cols ; ++j){
      activated(i, j) = max(0.0, matrix(i, j));
    }
  }
  return activated;
}

// Calculate the backwards pass for the activation.
// Parameters:
//  const Matrix& out: the activated output of the current layer.
//  const Matrix& prev_grad: the gradient from the next layer (towards the Loss).
// Returns:
//  Matrix: the gradients of this layer (to be passed to the previous layer).
Matrix backward_relu(const Matrix &out, const Matrix &prev_grad) {
  assert_same_size(prev_grad, out);
  Matrix grad = prev_grad;
  // TODO: Implement activation backward pass.
  for(int i = 0 ; i < grad.rows ; ++i){
    for(int j = 0 ; j < grad.cols ; ++j){
      grad(i, j) = (out(i, j) < 0) ? 0 : prev_grad(i, j);
    }
  }
  return grad;
}

// Calculate a Leaky ReLU activation.
// Use slope = 0.01
// Parameters:
//  Matrix &matrix: the input non-activated output of the layer.
// Returns:
Matrix forward_lrelu(const Matrix &matrix) {
  Matrix activated = matrix;
  // TODO: Implement forward activation.
  for(int i = 0 ; i < matrix.rows ; ++i){
    for(int j = 0 ; j < matrix.cols ; ++j){
      activated(i, j) = (matrix(i, j) < 0) ? 0.01 * matrix(i, j) : matrix(i, j);
    }
  }
  return activated;
}

// Calculate the backwards pass for the activation.
// Parameters:
//  const Matrix& out: the activated output of the current layer.
//  const Matrix& prev_grad: the gradient from the next layer (towards the Loss).
// Returns:
//  Matrix: the gradients of this layer (to be passed to the previous layer).
Matrix backward_lrelu(const Matrix &out, const Matrix &prev_grad) {
  assert_same_size(prev_grad, out);
  Matrix grad = prev_grad;
  // TODO: Implement activation backward pass.
  for(int i = 0 ; i < grad.rows ; ++i){
    for(int j = 0 ; j < grad.cols ; ++j){
      grad(i, j) = (out(i, j) < 0) ? 0.01 * prev_grad(i, j) : prev_grad(i, j);
    }
  }
  return grad;
}

// Calculate a Softmax activation.
// https://mattpetersen.github.io/softmax-with-cross-entropy
// Parameters:
//  Matrix &matrix: the input non-activated output of the layer.
// Returns:
Matrix forward_softmax(const Matrix &matrix) {
  Matrix activated = matrix;
  // TODO: Implement forward activation.
  for(int i = 0 ; i < matrix.rows ; ++i){

    const double* max_x = std::max_element(matrix[i], matrix[i] + matrix.cols);

    float sum_x = 0 ;
    for(int j = 0 ; j < matrix.cols ; ++j){
      activated(i, j) = std::exp(matrix(i, j) - *max_x);
      sum_x += activated(i, j);
    }

    for(int j = 0 ; j < matrix.cols ; ++j){
      activated(i, j) /= sum_x;
    }

  }
  return activated;
}

// Computes the Jacobian of the softmax function for a single row.
//
// Parameters:
//  Matrix &out_row: a 1xM vector matrix representing the output activation of a softmax function.
// Returns:
//  an MxM matrix representing the Jacobian matrix.
Matrix softmax_jacobian(const Matrix &out_row) {
  assert(out_row.rows == 1);
  Matrix jacobian(out_row.cols, out_row.cols);
  // TODO: Implement the Jacobian matrix.
  // Do whatever you want here, but here's some structure to get you started.
  for (int j = 0; j < out_row.cols; j++) {
    for (int k = 0; k < out_row.cols; k++) {
      jacobian(j, k) = (j == k) ? out_row(0, j) - out_row(0, j) * out_row(0, k): -out_row(0, j) * out_row(0, k);
    }
  }
  assert(jacobian.rows == out_row.cols);
  assert(jacobian.cols == out_row.cols);
  return jacobian;
}

// Computes the backwards pass for the softmax function.
Matrix backward_softmax(const Matrix &out, const Matrix &prev_grad) {
  assert_same_size(prev_grad, out);
  // TODO: Implement activation backward pass.
  Matrix grad = prev_grad;
  // Multiply previous gradient with Jacobian.
  for (int i = 0 ; i < out.rows ; ++i) {
    Matrix jacobian = softmax_jacobian(out.get_row(i));
    Matrix row_grad = prev_grad.get_row(i);
    // TODO: Implement the softmax backward pass.
    // cout << row_grad.rows << ' ' << row_grad.cols << endl;
    // cout << jacobian.rows << ' ' << jacobian.cols << endl;
    for(int j = 0 ; j < out.cols ; ++j){
      grad(i, j) = 0;
      for(int k = 0 ; k < jacobian.cols ; ++k){
        grad(i, j) += row_grad(0, k) * jacobian(j, k);
      }
    }
  }
  return grad;
}

// Run an activation function on each element in a matrix,
//
// Matrix& m: Input to activation function
// Activation a: function to run
// return the activated matrix
Matrix forward_activate_matrix(const Matrix &matrix, Activation a) {
  if (a == LINEAR) {
    return forward_linear(matrix);
  } else if (a == LOGISTIC) {
    return forward_logistic(matrix);
  } else if (a == TANH) {
    return forward_tanh(matrix);
  } else if (a == RELU) {
    return forward_relu(matrix);
  } else if (a == LRELU) {
    return forward_lrelu(matrix);
  } else if (a == SOFTMAX) {
    return forward_softmax(matrix);
  } else {
    assert(false); // Invalid activation.
  }
}

// Calculates the gradient of an activation function
// and multiplies it into the initial gradient for a layer
//
// const Matrix& out: an activated layer output
// Activation a: activation function for a layer
// Matrix& grad: before activation gradient (initial layer gradient)
// returns: Matrix that is after applying the activation gradien
Matrix backward_activate_matrix(const Matrix &out, const Matrix &grad, Activation a) {
  if (a == LINEAR) {
    return backward_linear(out, grad);
  } else if (a == LOGISTIC) {
    return backward_logistic(out, grad);
  } else if (a == TANH) {
    return backward_tanh(out, grad);
  } else if (a == RELU) {
    return backward_relu(out, grad);
  } else if (a == LRELU) {
    return backward_lrelu(out, grad);
  } else if (a == SOFTMAX) {
    return backward_softmax(out, grad);
  } else {
    assert(false); // Invalid activation.
  }
}
