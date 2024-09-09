#include <Rcpp.h>
#include <RcppEigen.h>
#include <memory>
#include <cmath>
#include <Eigen/Dense>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace Eigen;
using namespace RcppParallel;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins("cpp11")]]

// An abstract base class In_LFun
class In_LFun {
public:
  virtual double operator()(double eta) const = 0;
};

// Inverse link function for the Bernoulli distribution
class B_In_LFun : public In_LFun {
public:
  double operator()(double eta) const {
    return 1.0 / (1.0 + exp(-eta));
  }
};


// Inverse link function for the Poisson distribution
class P_In_LFun : public In_LFun {
public:
  double operator()(double eta) const {
    return exp(eta);
  }
};

// Base class for the W matrix
class Ca_W {
public:
  virtual Eigen::VectorXd operator()(Eigen::VectorXd mu) const = 0;
};

// W matrix for the Gaussian distribution
class G_Ca_W : public Ca_W {
public:
  Eigen::VectorXd operator()(Eigen::VectorXd mu) const {
    int n = mu.size();
    double v = (mu.array() - mu.mean()).square().sum() / (n - 1);
    double v1 = 1 / v;
    Eigen::VectorXd wii = Eigen::VectorXd::Constant(n, v1);
    return wii;
  }
};

// W matrix for the Bernoulli distribution
class B_Ca_W : public Ca_W {
public:
  Eigen::VectorXd operator()(Eigen::VectorXd mu) const {
    return mu.array() * (1 - mu.array());
  }
};

// W matrix for another Bernoulli distribution
class B2_Ca_W : public Ca_W {
public:
  Eigen::VectorXd nn;

  // Constructor accepting nn vector as parameter and initializing
  B2_Ca_W(const Eigen::VectorXd& nn) : nn(nn) {}

  Eigen::VectorXd operator()(Eigen::VectorXd mu) const {
    int n = mu.size();
    Eigen::VectorXd one_minus_mu = Eigen::VectorXd::Constant(n, 1.0) - mu.cwiseQuotient(nn);
    Eigen::VectorXd wii = mu.cwiseProduct(one_minus_mu);
    return wii;
  }
};

// W matrix for the Poisson distribution
class P_Ca_W : public Ca_W {
public:
  Eigen::VectorXd operator()(Eigen::VectorXd mu) const {
    return mu;
  }
};

// W matrix for the Negative binomial distribution
class N_Ca_W : public Ca_W {
public:
  double a;
  Eigen::VectorXd y;

  N_Ca_W(double a) :  a(a) {}

  Eigen::VectorXd operator()(Eigen::VectorXd mu) const {

    Eigen::VectorXd W = mu.array() / (1 + a * mu.array());

    return W;
  }
};

// Base class for the A matrix
class Ca_A {
public:
  virtual Eigen::VectorXd operator()(Eigen::VectorXd y, Eigen::VectorXd mu) const = 0;
};

// A matrix for the distribution
class B_Ca_A : public Ca_A {
public:
  Eigen::VectorXd operator()(Eigen::VectorXd y, Eigen::VectorXd mu) const {
    return y - mu;
  }
};

// A matrix for the Negative binomial distribution
class N_Ca_A : public Ca_A {
public:
  double a;

  N_Ca_A(double a) : a(a) {}

  Eigen::VectorXd operator()(Eigen::VectorXd y, Eigen::VectorXd mu) const {
    return (y - mu).array() / (1 + a * mu.array());
  }
};

//  Summary:
//   Parameters estimation using the Newton-Raphson method
//
//  Parameters:
//
//   y: response variable vector
//   nn: weight vector
//   X: covariate matrix.
//   family: the distribution family
//   a: theta parameter
//
//  Return:
//   estimated parameter values (beta)
//
// [[Rcpp::export]]
Eigen::VectorXd glm_NR(Eigen::VectorXd y,Eigen::VectorXd nn, Eigen::MatrixXd X,
                       std::string family,double a) {
  int maxIter = 500;
  double tol = 1e-6;
  int n = X.rows();
  int p = X.cols();
  int iter;
  // Eigen::MatrixXd nX(n, p + 1);
  // nX << Eigen::MatrixXd::Ones(n, 1), X;

  Eigen::VectorXd beta(p);
  beta.setZero();
  if (family == "gaussian") {
    Eigen::MatrixXd XtX = X.transpose() * X;
    Eigen::VectorXd Xty = X.transpose() * y;

    Eigen::VectorXd beta = XtX.fullPivLu().solve(Xty);
    return beta;

  } else if (family == "binomial") {

    if (!nn.isZero() && y.maxCoeff() > 1) {
      double mean_y = y.sum() / nn.sum();
      beta[0] = std::log(mean_y / (1 - mean_y));
    } else {
      double meanY = y.mean();
      beta[0] = log(meanY / (1 - meanY));
    }

  } else if (family == "poisson" || family == "negative.binomial") {
    beta[0] = log(y.mean());
  } else {
    Rcpp::stop("Unsupported family type");
  }
  //std::cout << "beta: " << beta.transpose() << std::endl;

  // Select the inverse connection function and formula calculation function based on family
  std::unique_ptr<In_LFun> inlfunc = nullptr;
  std::unique_ptr<Ca_W> cawc = nullptr;
  std::unique_ptr<Ca_A> caac = nullptr;

  if(family == "binomial") {
    inlfunc = std::unique_ptr<B_In_LFun>(new B_In_LFun());
    caac = std::unique_ptr<B_Ca_A>(new B_Ca_A());

    if (!nn.isZero() && y.maxCoeff() > 1) {
      cawc = std::unique_ptr<B2_Ca_W>(new B2_Ca_W(nn));
    } else {
      cawc = std::unique_ptr<B_Ca_W>(new B_Ca_W());
    }

  } else if (family == "poisson") {
    inlfunc = std::unique_ptr<P_In_LFun>(new P_In_LFun());
    cawc = std::unique_ptr<P_Ca_W>(new P_Ca_W());
    caac = std::unique_ptr<B_Ca_A>(new B_Ca_A());
  } else if (family == "negative.binomial") {
    inlfunc = std::unique_ptr<P_In_LFun>(new P_In_LFun());
    cawc = std::unique_ptr<N_Ca_W>(new N_Ca_W(a));
    caac = std::unique_ptr<N_Ca_A>(new N_Ca_A(a));
  } else {
    Rcpp::stop("Unsupported family type");
  }
  // Execute iteration process
  for (int iter = 0; iter < maxIter; ++iter) {
    // Calculate eta and mu vectors based on the link function and inverse link function

    Eigen::VectorXd eta = X * beta;
    Eigen::VectorXd mu(n);

    if (family == "binomial" && (!nn.isZero() && y.maxCoeff() > 1)) {
      for (int i = 0; i < n; ++i) {
        mu[i] = (*inlfunc)(eta[i]) * nn[i];
      }
    }else{
      for (int i = 0; i < n; ++i) {
        mu[i] = (*inlfunc)(eta[i]);
      }
    }

    Eigen::VectorXd W = (*cawc)(mu);

    Eigen::VectorXd A = (*caac)(y, mu);

    Eigen::MatrixXd J = X.transpose() * W.asDiagonal() * X;
    //std::cout << "j:\n " << J << std::endl;
    Eigen::VectorXd U = X.transpose() * A;
    // Update the parameter vector beta
    Eigen::VectorXd nbeta = beta + J.fullPivLu().solve(U);
    //std::cout << "nbeta: " << nbeta.transpose() << std::endl;
    // Check for convergence
    //double diff = sqrt(dot(nbeta - beta,nbeta - beta) / dot(nbeta,beta));
    double diff = (nbeta - beta).norm() / beta.norm();
    //std::cout << "diff " << diff << std::endl;
    if (diff < tol) {
      break;
    }

    // Update the parameter vector beta
    beta = nbeta;
  }


  // Check if the iteration did not converge
  if (iter == maxIter) {
    Rcpp::warning("Newton-Raphson iteration did not converge.");
  }

  return beta;
}

// Summary:
//  The expected mu vector
//
// Parameters:
//   y: response variable vector
//   nn: weight vector
//   X: covariate matrix
//   family: the distribution family
//   beta: parameter estimate vector
//
// Return:
//   expected mu vector
//
// [[Rcpp::export]]
Eigen::VectorXd Ca_Mu(Eigen::VectorXd  y,Eigen::VectorXd nn, Eigen::MatrixXd X, std::string family,
                      Eigen::VectorXd  beta) {
  int n = X.rows();
  int p = X.cols();

  // Add a column of all 1’s to the first column of the X matrix
  // NumericMatrix newX(n, p + 1);
  // newX(_, 0) = 1.0;
  // for (int j = 0; j < p; ++j) {
  //   newX(_, j + 1) = X(_, j);
  // }
  // Eigen::MatrixXd nX(n, p + 1);
  // nX.col(0).setOnes();
  // nX.block(0, 1, n, p) = X;

  // Create an inverse join function object
  std::unique_ptr<In_LFun> inlfunc = nullptr;
  if (family == "gaussian") {
    //inlfunc = new G_In_LFun();
    Eigen::VectorXd mu = X * beta;
    return mu;

  } else if (family == "binomial") {
    inlfunc = std::unique_ptr<B_In_LFun>(new B_In_LFun());
  } else if (family == "poisson" || family == "negative.binomial") {
    inlfunc = std::unique_ptr<P_In_LFun>(new P_In_LFun());
  } else {
    throw std::runtime_error("Invalid family type");
  }

  // Calculate mu vector
  Eigen::VectorXd eta = X * beta;
  Eigen::VectorXd mu(n);

  for (int i = 0; i < n; ++i) {
    mu[i] = (*inlfunc)(eta[i]);
  }

  return mu;
}

// When family is negative binomial distribution, calculate the dispersion parameter
double pcss(Eigen::VectorXd y, Eigen::VectorXd mu, int p) {
  int N = y.size();
  double fi = 0.0;

  for (int i = 0; i < N; ++i) {
    double nr = pow(y[i] - mu[i], 2);
    double dr = mu[i] * (N - p);
    fi += nr / dr;
  }

  return fi;
}

// Summary:
// The diagonal elements of the W matrix
//
// Parameters:
//   y: response variable vector
//   nn: weight vector
//   family: the distribution family
//   mu: expected vector
//   beta: estimated parameter values
//   a: theta parameter
//
// Returns:
//   Vector of diagonal elements of the W matrix (Wii)
//
// [[Rcpp::export]]
Eigen::VectorXd Ca_Wii(Eigen::VectorXd y,Eigen::VectorXd nn, std::string family,
                       Eigen::VectorXd mu,Eigen::VectorXd beta,double a) {
  int n = y.size();
  int p =beta.size();

  // Create W matrix calculation function object
  std::unique_ptr<Ca_W> cawc = nullptr;
  if (family == "gaussian") {
    // int n = y.size();
    // double v = (y.array() - y.mean()).square().sum() / (n - 1);
    // double v1 = 1 / v;
    // Eigen::VectorXd Wii = Eigen::VectorXd::Constant(n, v1);
    //cawc = std::unique_ptr<G_Ca_W>(new G_Ca_W());
    Eigen::VectorXd Wii = Eigen::VectorXd::Ones(n);
    return Wii;

  } else if (family == "binomial") {
    cawc = std::unique_ptr<B_Ca_W>(new B_Ca_W());
  } else if (family == "poisson") {
    cawc = std::unique_ptr<P_Ca_W>(new P_Ca_W());
  } else if (family == "negative.binomial") {
    cawc = std::unique_ptr<N_Ca_W>(new N_Ca_W(a));
  }else {
    throw std::runtime_error("Invalid family type");
  }

  // Calculate W diagonal element vector
  Eigen::VectorXd Wii = (*cawc)(mu);
  // if (family == "poisson") {
  //   double fi = pcss(y, mu, p);
  //   Wii /= fi;
  // }
  return Wii;
}

// Summary:
// Generate a submatrix by replicating the matrix based on envir
//
// Parameters:
//   X: input matrix
//   envir: number of environmental factors
//
// Returns:
//   Transformed matrix with replicated blocks (Xs)
//
// [[Rcpp::export]]
Eigen::MatrixXd generateSubmatrix(const Eigen::MatrixXd& X, int envir) {
  int n = X.rows();
  int m = X.cols();
  Eigen::MatrixXd Xs(n * envir, m * envir);
  Xs.setZero();
  for (int j = 0; j < envir; j++) {
    Xs.block(j * n, j * m, n, m) = X;
  }

  return Xs;
}

// Model testing, score testing
double score(Eigen::VectorXd y,Eigen::VectorXd nn, Eigen::MatrixXd X, std::string family,
             Eigen::VectorXd beta,double a) {
  int n = X.rows();
  int p = X.cols();

  // Eigen::MatrixXd nX(n, p + 1);
  // nX.col(0).setOnes();
  // nX.block(0, 1, n, p) = X;

  Eigen::VectorXd mu = Ca_Mu(y,nn, X, family, beta);
  Eigen::VectorXd Wii = Ca_Wii(y,nn, family, mu,beta,a);
  std::unique_ptr<Ca_A> caac = nullptr;
  if(family == "negative.binomial") {
    caac = std::unique_ptr<N_Ca_A>(new N_Ca_A(a));
  } else {
    caac = std::unique_ptr<B_Ca_A>(new B_Ca_A());
  }
  Eigen::VectorXd A = (*caac)(y, mu);
  // if (family == "poisson") {
  //   double fi = pcss(y, mu, p);
  //   A /= fi;
  // }

  // Calculate the score statistic U and the information matrix J
  Eigen::VectorXd U = X.transpose() * A;
  Eigen::MatrixXd J = X.transpose() * Wii.asDiagonal() * X;

  // Calculate the score values
  double score = U.transpose() * J.inverse() * U;

  if (family == "poisson") {
    double fi = pcss(y, mu, p);
    score /= fi;
  }

  return score;
}

// Standard Error
Eigen::VectorXd SE(Eigen::VectorXd y,Eigen::VectorXd nn, Eigen::MatrixXd X, std::string family,
                   Eigen::VectorXd beta,double a) {
  int n = X.rows();
  int p = X.cols();

  // Eigen::MatrixXd nX(n, p + 1);
  // nX.col(0).setOnes();
  // nX.block(0, 1, n, p) = X;

  Eigen::VectorXd mu = Ca_Mu(y,nn, X, family, beta);
  Eigen::VectorXd Wii = Ca_Wii(y,nn, family, mu,beta,a);

  Eigen::MatrixXd J = X.transpose() * Wii.asDiagonal() * X;

  Eigen::VectorXd se = J.inverse().diagonal().cwiseSqrt();

  return se;
}

//// Model testing, wald testing
Eigen::VectorXd wald(Eigen::VectorXd beta, Eigen::VectorXd se) {

  Eigen::VectorXd w = (beta.array() / se.array()).square();

  return w;
}

// Summary:
// Forward selection method
//
// Parameters:
//   y: response variable vector
//   nn: weight vector
//   X: covariate matrix
//   G: genotype matrix
//   family: the distribution family
//   envir: number of environmental factors
//   threshold: selection threshold
//   a: theta parameter.
//   x_indices: indices of features to be selected
//   y_indices: indices of features to be retained
//
// Returns:
//   Vector of indices of selected significant features (x_selected)
//
// [[Rcpp::export]]
Eigen::VectorXi forward_selection(const Eigen::VectorXd& y,const Eigen::VectorXd& nn,
                                  const Eigen::MatrixXd& X,const Eigen::MatrixXd& G,
                                  const std::string& family, int envir,
                                  double threshold, double a , Eigen::VectorXi& x_indices,
                                  Eigen::VectorXi& y_indices) {
  int n = G.rows();
  int m = G.cols();
  int c = X.cols();
  int d = x_indices.size();
  int b = y_indices.size();
  //double threshold = 0.5 / m;
  Eigen::VectorXi& x_selected = x_indices;
  Eigen::VectorXi& y_remaining = y_indices;
  if (x_indices == y_indices) {
    return x_selected;
  } else {

    Eigen::MatrixXd nX(n * envir, d * envir + c);

    Eigen::MatrixXd Xs(n * envir, d * envir);
    Eigen::MatrixXd G_selected(n, d);
    for (int j = 0; j < d; j++) {
      G_selected.col(j) = G.col(x_selected(j));
    }

    Xs = generateSubmatrix(G_selected, envir);

    nX << X, Xs;

    while (b > 0 ) {
      //int a = x_selected.size();
      double min_p_value = std::numeric_limits<double>::max();
      int min_p_index = -1;

      Eigen::VectorXd beta = glm_NR(y,nn, nX, family,a);
      beta.conservativeResize(beta.size() + envir);
      beta.tail(envir).setZero();

      for (int i = 0; i < b; i++) {
        int index = y_remaining(i);
        Eigen::VectorXd G_col = G.col(index);
        Eigen::MatrixXd Xss = generateSubmatrix(G_col, envir);
        Eigen::MatrixXd Xcs(n * envir, d * envir + c + envir);
        Xcs << nX, Xss;
        double s = score(y,nn, Xcs, family, beta,a);

        double p_value = R::pchisq(s, envir, false, false);

        if (p_value < min_p_value) {
          min_p_value = p_value;
          min_p_index = i;
        }
      }

      if (min_p_value <= threshold) {
        int selected_index = y_remaining(min_p_index);
        x_selected.conservativeResize(d + 1);
        x_selected(d) = selected_index;

        y_remaining.segment(min_p_index, b - min_p_index - 1) = y_remaining.segment(min_p_index + 1, b - min_p_index - 1);
        y_remaining.conservativeResize(b - 1);

        nX.conservativeResize(n * envir, d * envir + c + envir);
        nX.block(0, d * envir + c, n * envir, envir) = generateSubmatrix(G.col(selected_index), envir);

        d++;
        b--;
      } else {
        break;
      }
    }
    return x_selected;
  }
}

// When performing a two-dimensional scan, generate an index into the genotype matrix
std::vector<std::tuple<int, int>> getIndices(int m) {
  std::vector<std::tuple<int, int>> indices;
  indices.reserve(m * (m - 1) / 2);

  for (int i = 0; i < m - 1; i++) {
    for (int j = i + 1; j < m; j++) {
      indices.push_back(std::make_tuple(i, j));
    }
  }

  return indices;
}

// Custom worker class for parallel computation
struct ScoreTestWorker : public Worker {
  const Eigen::VectorXd y;
  const Eigen::VectorXd mu;
  const Eigen::VectorXd Wii;
  const Eigen::MatrixXd G;
  const Eigen::VectorXd u;
  const Eigen::MatrixXd& proj;
  const double b;
  const int envir;
  bool epistasis_test;
  Eigen::VectorXd& x2;
  Eigen::VectorXd& p_value;

  //std::size_t begin;
  //std::size_t end;
  ScoreTestWorker(const Eigen::VectorXd& y, const Eigen::VectorXd& mu,
                  const Eigen::VectorXd& Wii,const Eigen::MatrixXd& G,
                  const Eigen::VectorXd& u,const Eigen::MatrixXd& proj,
                  const double b, const int envir,bool epistasis_test,
                  Eigen::VectorXd& x2, Eigen::VectorXd& p_value)
    : y(y), mu(mu), Wii(Wii), G(G),u(u),proj(proj),b(b), envir(envir), epistasis_test(epistasis_test),
      x2(x2), p_value(p_value) {}

  void operator()(std::size_t begin, std::size_t end) {
    // int numIterations = end - begin;
    // lor.resize(numIterations, envir);
    // odds_ratio.resize(numIterations, envir);
    // se.resize(numIterations, envir);
    // z_value.resize(numIterations, envir);

    if (!epistasis_test) {

      for (std::size_t i = begin; i < end; i++) {
        VectorXd X = G.col(i);
        MatrixXd Xs = generateSubmatrix(X,envir);
        //MatrixXd Xs = MatrixXd::Zero(y.size(), envir); // Xs变为(envir*n)×envir矩阵
        // for (int j = 0; j < envir; j++) {
        //   Xs.col(j) = G.col(i).segment(j * y.size(), y.size()); // Xs的每一列为对应环境的G列
        // }
        //
        VectorXd Tscore = Xs.transpose() * u; // Tscore向量是envir×1
        MatrixXd varT = Xs.transpose() * Wii.asDiagonal() * Xs; // varT 是envir×envir
        MatrixXd invVarT = varT.inverse();
        x2(i) = Tscore.transpose() * invVarT * Tscore;
        p_value(i) = R::pchisq(x2(i), envir, false, false);
        if(p_value(i) <= b){
          Eigen::MatrixXd Xs_proj = Xs - proj * Xs ;
          Eigen::VectorXd Tscore = Xs_proj.transpose() * u;
          Eigen::MatrixXd varT = Xs_proj.transpose() * Wii.asDiagonal() * Xs_proj;
          Eigen::MatrixXd invVarT = varT.inverse();
          x2(i) = Tscore.transpose() * invVarT * Tscore;
          p_value(i) = R::pchisq(x2(i), envir, false, false);
        }

      }

    }else{
      int m = G.cols();
      std::vector<std::tuple<int, int>> indices = getIndices(m);

      for (std::size_t idx = begin; idx < end ; idx++) {
        int i = std::get<0>(indices[idx]);
        int j = std::get<1>(indices[idx]);
        Eigen::VectorXd X = G.col(i).array() * G.col(j).array();

        MatrixXd Xs = generateSubmatrix(X, envir);
        VectorXd Tscore = Xs.transpose() * u;
        MatrixXd varT = Xs.transpose() * Wii.asDiagonal() * Xs;
        MatrixXd invVarT = varT.inverse();
        x2(idx) = Tscore.transpose() * invVarT * Tscore;
        p_value(idx) = R::pchisq(x2(idx), envir, false, false);
        if(p_value(idx) <= b){
          Eigen::MatrixXd Xs_proj = Xs - proj * Xs ;
          Eigen::VectorXd Tscore = Xs_proj.transpose() * u;
          Eigen::MatrixXd varT = Xs_proj.transpose() * Wii.asDiagonal() * Xs_proj;
          Eigen::MatrixXd invVarT = varT.inverse();
          x2(idx) = Tscore.transpose() * invVarT * Tscore;
          p_value(idx) = R::pchisq(x2(idx), envir, false, false);
        }

      }

    }
  }
};

// Summary:
// Score test for model testing
//
// Parameters:
//   y: response variable vector
//   mu: expected vector
//   Wii: vector of diagonal elements of the W matrix
//   G: genotype matrix
//   proj: projection matrix
//   beta: estimated parameter values
//   family: the distribution family
//   envir: number of environmental factors
//   b: threshold for computation
//   a: theta parameter
//   epistasis_test: flag for conducting 2D scan
//   parallel: flag for parallel computation
//
// Returns:
// List containing 'x2' (score value) and 'p_value' (p-value)
//
// [[Rcpp::export]]
List score_test(Eigen::VectorXd y, Eigen::VectorXd mu, Eigen::VectorXd Wii,
                Eigen::MatrixXd G,Eigen::MatrixXd proj,Eigen::VectorXd beta,
                const std::string family,int envir = 1,double b = 0.01, double a = 0.0,
                bool epistasis_test = false, bool parallel = false) {
  int m = G.cols();
  int n = y.size();
  Eigen::VectorXd x2(m);
  Eigen::VectorXd p_value(m);
  Eigen::VectorXd u(n);
  if(family == "negative.binomial") {
    u = (y - mu).array() / (1 + a * mu.array());
  } else if(family == "poisson"){
    int p = beta.size();
    u = y - mu;
    double fi = pcss(y, mu, p);
    u /= sqrt(fi);
  }else {
    u = y - mu;
  }
  if ((proj.array() == 0).all()) {
    b = 0.0;
  }

  if (epistasis_test) {
    int num = m * (m - 1) / 2;
    x2.resize(num);
    p_value.resize(num);

    if (!parallel ) {
      int interaction_index = 0;
      for (int i = 0; i < m - 1; i++) {
        for (int j = i + 1; j < m; j++) {
          Eigen::VectorXd X = G.col(i).array() * G.col(j).array();
          Eigen::MatrixXd Xs = generateSubmatrix(X, envir);
          Eigen::VectorXd Tscore = Xs.transpose() * u;
          Eigen::MatrixXd varT = Xs.transpose() * Wii.asDiagonal() * Xs;
          Eigen::MatrixXd invVarT = varT.inverse();
          x2(interaction_index) = Tscore.transpose() * invVarT * Tscore;
          p_value(interaction_index) = R::pchisq(x2(interaction_index), envir, false, false);
          if(p_value(interaction_index) <= b){
            Eigen::MatrixXd Xs_proj = Xs - proj * Xs ;
            Eigen::VectorXd Tscore = Xs_proj.transpose() * u;
            Eigen::MatrixXd varT = Xs_proj.transpose() * Wii.asDiagonal() * Xs_proj;
            Eigen::MatrixXd invVarT = varT.inverse();
            x2(interaction_index) = Tscore.transpose() * invVarT * Tscore;
            p_value(interaction_index) = R::pchisq(x2(interaction_index), envir, false, false);
          }
          interaction_index++;
        }
      }

    } else {
      ScoreTestWorker worker(y, mu, Wii, G,u,proj,b,envir,true, x2, p_value);
      parallelFor(0, num, worker);
    }
  } else {
    if (!parallel) {
      for (int i = 0; i < m; i++) {
        Eigen::VectorXd X = G.col(i);
        Eigen::MatrixXd Xs = generateSubmatrix(X, envir);
        Eigen::VectorXd Tscore = Xs.transpose() * u;
        Eigen::MatrixXd varT = Xs.transpose() * Wii.asDiagonal() * Xs;
        Eigen::MatrixXd invVarT = varT.inverse();
        x2(i) = Tscore.transpose() * invVarT * Tscore;
        p_value(i) = R::pchisq(x2(i), envir, false, false);
        if(p_value(i) <= b){
          Eigen::MatrixXd Xs_proj = Xs - proj * Xs ;
          Eigen::VectorXd Tscore = Xs_proj.transpose() * u;
          Eigen::MatrixXd varT = Xs_proj.transpose() * Wii.asDiagonal() * Xs_proj;
          Eigen::MatrixXd invVarT = varT.inverse();
          x2(i) = Tscore.transpose() * invVarT * Tscore;
          p_value(i) = R::pchisq(x2(i), envir, false, false);
        }
      }
    } else{
      ScoreTestWorker worker(y, mu, Wii, G, u,proj,b,envir,false,x2, p_value);
      parallelFor(0, m, worker);
    }
  }

  return List::create(Named("x2") = x2,
                      Named("p_value") = p_value);

}

// RCPP_MODULE(glmModule) {
//   function("generateSubmatrix", &generateSubmatrix);
//   function("score_test", &score_test);
// }
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
// Export main function and extra functions in Rcpp module
// RCPP_MODULE(glmModule) {
//   function("glm_NR", &glm_NR);
//   function("Ca_Mu", &Ca_Mu);
//   function("Ca_Wii", &Ca_Wii);
// }
