#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#include <queue>

#include "optimalityfunctions.h"
#include "nullify_alg.h"

using namespace Rcpp;

// [[Rcpp::export]]
List genOptimalDesignCoord(Eigen::MatrixXd initialdesign, const Eigen::MatrixXd& candidatelist,
                          const std::string condition,
                          const Eigen::MatrixXd& momentsmatrix, Eigen::VectorXd initialRows,
                          Eigen::MatrixXd aliasdesign,
                          const Eigen::MatrixXd& aliascandidatelist,
                          double minDopt, double tolerance, int augmentedrows, int kexchange) {
  RNGScope rngScope;
  int nTrials = initialdesign.rows();
  double numberrows = initialdesign.rows();
  double numbercols = initialdesign.cols();
  int maxSingularityChecks = nTrials*100;
  int totalPoints = candidatelist.rows();
  Eigen::VectorXd candidateRow(nTrials);
  candidateRow.setZero();
  Eigen::MatrixXd test(initialdesign.cols(), initialdesign.cols());
  test.setZero();
  if(nTrials < candidatelist.cols()) {
    throw std::runtime_error("Too few runs to generate initial non-singular matrix: increase the number of runs or decrease the number of parameters in the matrix");
  }
  //Check for singularity from a column perfectly correlating with the intercept.
  for(int j = 1; j < candidatelist.cols(); j++) {
    if(candidatelist.col(0).cwiseEqual(candidatelist.col(j)).all()) {
      throw std::runtime_error("Singular model matrix from factor aliased into intercept, revise model");
    }
  }
  Eigen::VectorXi shuffledindices;
  //Checks if the initial matrix is singular. If so, randomly generates a new design maxSingularityChecks times.
  for (int check = 0; check < maxSingularityChecks; check++) {
    if(!isSingular(initialdesign)) {
      break; //design is nonsingular
    }
    if(nTrials <= totalPoints) {
      shuffledindices = sample_noreplace(totalPoints, nTrials);
    } else {
      shuffledindices = sample_noreplace(nTrials, nTrials);
      for(int i = 0; i < shuffledindices.size(); i++) {
        shuffledindices(i) %= totalPoints;
      }
    }

    for (int i = augmentedrows; i < nTrials; i++) {
      initialdesign.row(i) = candidatelist.row(shuffledindices(i));
      aliasdesign.row(i) = aliascandidatelist.row(shuffledindices(i));
      initialRows(i) = shuffledindices(i) + 1; //R indexes start at 1
    }
  }
  //If initialdesign is still singular, use the Gram-Schmidt orthogonalization procedure, which
  //should return a non-singular matrix if one can be constructed from the candidate set
  if (isSingular(initialdesign)) {
    Eigen::VectorXi initrows = orthogonal_initial(candidatelist, nTrials);
    //If all elements are equal here, nullification algorithm was unable to find a design--return NA
    if(initrows.minCoeff() == initrows.maxCoeff()) {
      return(List::create(_["indices"] = NumericVector::get_na(), _["modelmatrix"] = NumericMatrix::get_na(), _["criterion"] = NumericVector::get_na()));
    }

    //Replace non-augmented rows with orthogonal design
    for (int i = 0; i < nTrials - augmentedrows; i++) {
      initialdesign.row(i + augmentedrows) = candidatelist.row(initrows(i));
      aliasdesign.row(i + augmentedrows) = aliascandidatelist.row(initrows(i));
      initialRows(i + augmentedrows) = initrows(i) + 1; //R indexes start at 1
    }

    //Shuffle design
    Eigen::VectorXi initrows_shuffled = sample_noreplace(nTrials - augmentedrows, nTrials - augmentedrows);
    for (int i = augmentedrows; i < nTrials; i++) {
      initialdesign.row(i) = initialdesign.row(augmentedrows + initrows_shuffled(i));
      aliasdesign.row(i) = aliasdesign.row(augmentedrows + initrows_shuffled(i));
      initialRows(i) = augmentedrows + initrows_shuffled(i) + 1; //R indexes start at 1
    }
  }
  //If still no non-singular design, returns NA.
  if (isSingular(initialdesign)) {
    return(List::create(_["indices"] = NumericVector::get_na(), _["modelmatrix"] = NumericMatrix::get_na(), _["criterion"] = NumericVector::get_na()));
  }

  bool found = true;
  double del = 0;
  int entryx = 0;
  int entryy = 0;
  double newOptimum = 0;
  double priorOptimum = 0;
  double minDelta = tolerance;
  double newdel;
  double xVx;

  //Initialize matrices for rank-2 updates.
  Eigen::MatrixXd identitymat(2,2);
  identitymat.setIdentity(2,2);
  Eigen::MatrixXd f1(initialdesign.cols(),2);
  Eigen::MatrixXd f2(initialdesign.cols(),2);
  Eigen::MatrixXd f2vinv(2,initialdesign.cols());

  //Transpose matrices for faster element access
  Eigen::MatrixXd initialdesign_trans = initialdesign.transpose();
  Eigen::MatrixXd candidatelist_trans = candidatelist.transpose();
  Eigen::MatrixXd V = (initialdesign.transpose()*initialdesign).partialPivLu().inverse();
  //Generate a D-optimal design
  if(condition == "D") {
    newOptimum = 1.0;
    priorOptimum = newOptimum/2;

    while((newOptimum - priorOptimum)/priorOptimum > minDelta) {
      priorOptimum = newOptimum;
      //Calculate k-exchange coordinates
      std::priority_queue<std::pair<double, int>> q;
      float min_val = -INFINITY;
      int k = kexchange - augmentedrows;
      if(kexchange != nTrials) {
        for (int i = augmentedrows; i < nTrials; i++) {
          float temp_val = -initialdesign_trans.col(i).transpose() * V * initialdesign_trans.col(i);
          if(temp_val == min_val) {
            k++;
          } else if(temp_val > min_val) {
            min_val = temp_val;
            k = kexchange - augmentedrows;
          }
          q.push(std::pair<double, int>(temp_val, i));
        }
      } else {
        for (int i = augmentedrows; i < nTrials; i++) {
          q.push(std::pair<double, int>(-i, i));
        }
      }

      for (int j = 0; j < k; j++) {
        Rcpp::checkUserInterrupt();
        int i = q.top().second;
        q.pop();
        found = false;
        entryy = 0;
        del=0;
        xVx = initialdesign_trans.col(i).transpose() * V * initialdesign_trans.col(i);
        //Search through all candidate set points to find best switch (if one exists).

        search_candidate_set(V, candidatelist_trans, initialdesign_trans.col(i), xVx, entryy, found, del);

        if (found) {
          //Update the inverse with the rank-2 update formula.
          rankUpdate(V,initialdesign_trans.col(i),candidatelist_trans.col(entryy),identitymat,f1,f2,f2vinv);

          //Exchange points and re-calculate current criterion value.
          initialdesign_trans.col(i) = candidatelist_trans.col(entryy);
          candidateRow[i] = entryy+1;
          initialRows[i] = entryy+1;
          newOptimum = newOptimum * (1 + del);
        } else {
          candidateRow[i] = initialRows[i];
        }
      }
    }
    initialdesign = initialdesign_trans.transpose();
    newOptimum = calculateDEff(initialdesign,numbercols,numberrows);
    if(std::isinf(newOptimum)) {
      newOptimum = calculateDEffLog(initialdesign,numbercols,numberrows);
    }
  }

  //return the model matrix and a list of the candidate list indices used to construct the run matrix
  return(List::create(_["indices"] = candidateRow, _["modelmatrix"] = initialdesign, _["criterion"] = newOptimum));
}

