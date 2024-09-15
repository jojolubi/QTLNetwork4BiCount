#ifndef PAR_GLM_H
#define PAR_GLM_H

#include <Rcpp.h>
#include <RcppEigen.h>
#include <memory>
#include <cmath>
#include <Eigen/Dense>
#include <RcppParallel.h>

using namespace Rcpp;
using namespace Eigen;
using namespace RcppParallel;

Eigen::VectorXd glm_NR(Eigen::VectorXd y, Eigen::VectorXd nn, Eigen::MatrixXd X,
                       std::string family, double a);

Eigen::VectorXd Ca_Mu(Eigen::VectorXd y, Eigen::VectorXd nn, Eigen::MatrixXd X,
                      std::string family, Eigen::VectorXd beta);

double pcss(Eigen::VectorXd y, Eigen::VectorXd mu, int p);

Eigen::VectorXd Ca_Wii(Eigen::VectorXd y, Eigen::VectorXd nn, std::string family,
                       Eigen::VectorXd mu, Eigen::VectorXd beta, double a);

Eigen::MatrixXd generateSubmatrix(const Eigen::MatrixXd& X, int envir);

double score(Eigen::VectorXd y,Eigen::VectorXd nn, Eigen::MatrixXd X, std::string family,
             Eigen::VectorXd beta,double a);

Eigen::VectorXd SE(Eigen::VectorXd y,Eigen::VectorXd nn, Eigen::MatrixXd X, std::string family,
                   Eigen::VectorXd beta,double a);

Eigen::VectorXd wald(Eigen::VectorXd beta, Eigen::VectorXd se);

Eigen::VectorXi forward_selection(const Eigen::VectorXd& y, const Eigen::VectorXd& nn,
                                  const Eigen::MatrixXd& X, const Eigen::MatrixXd& G,
                                  const std::string& family, int envir, double threshold,
                                  double a, Eigen::VectorXi& x_indices, Eigen::VectorXi& y_indices);

std::vector<std::tuple<int, int>> getIndices(int m);

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
  
  ScoreTestWorker(const Eigen::VectorXd& y, const Eigen::VectorXd& mu,
                  const Eigen::VectorXd& Wii, const Eigen::MatrixXd& G,
                  const Eigen::VectorXd& u, const Eigen::MatrixXd& proj,
                  const double b, const int envir, bool epistasis_test,
                  Eigen::VectorXd& x2, Eigen::VectorXd& p_value)
    : y(y), mu(mu), Wii(Wii), G(G), u(u), proj(proj),
      b(b), envir(envir), epistasis_test(epistasis_test),
      x2(x2), p_value(p_value) {}

  void operator()(std::size_t begin, std::size_t end);
  
};

List score_test(Eigen::VectorXd y, Eigen::VectorXd mu, Eigen::VectorXd Wii,
                Eigen::MatrixXd G, Eigen::MatrixXd proj, Eigen::VectorXd beta,
                const std::string family, int envir, double b, double a,
                bool epistasis_test = false, bool parallel = false);

#endif //PAR_GLM_H