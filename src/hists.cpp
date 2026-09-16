// bacon_hist.cpp - calculates ages and histograms for depths based on provided data
// does not read/write files - instead everything resides within the R session

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>

using namespace Rcpp;



// to calculate quantiles with interpolation as done in R
// (R offers 9 types of quantile calculations, with default type 7)
double quantileR(const NumericVector& x, double p) {
    int n = x.size();
    if (n == 0) return NA_REAL;

    // R type-7 quantile definition, m = 1-p, h = (n - 1) * p + 1
    double h = 1 + (n - 1) * p;
    int k = floor(h);
    double gm = h - k;

    // edge cases
    if (k <= 1) return x[0];
    if (k >= n) return x[n - 1];

    // linear interpolation
    return (1 - gm) * x[k - 1] + gm * x[k];
}



// [[Rcpp::export]]
NumericMatrix depths_ageranges(
  const NumericVector& depths,
  const NumericMatrix& out,
  const NumericVector& elbows, // elbow depths: length K+1
  int n_rows = 4000,
  double prob = 0.95) {

  int n_depths = depths.size(); // number of depths to calculate
  int K = elbows.size() - 1; // number of sections
  double thick = elbows[1] - elbows[0]; // length of sections
  
  NumericMatrix theta(n_rows, K + 1); // age
  NumericMatrix quantiles(n_depths, 5); // depths, min95%, max95%, median, mean
  double lower = (1.0 - prob) / 2.0; // quantile 0.025 if prob=0.95
  double upper = 1.0 - lower; // quantile 0.975 if prob=0.95

  for (int r = 0; r < n_rows; r++) { // find elbow ages
    theta(r, 0) = out(r, 0); // core-top age for iteration r
    for (int k = 1; k <= K; k++) 
      theta(r, k) = theta(r, k - 1) + out(r, k) * thick;
  }

  // precompute section index k for each depth
  std::vector<int> depth_section(n_depths);
  for (int i = 0; i < n_depths; i++) {
    double d = depths[i];
    int idx = std::upper_bound(elbows.begin(), elbows.end(), d)
      - elbows.begin() - 1;

    if (idx < 0) idx = 0;
    if (idx >= K) idx = K - 1;
    depth_section[i] = idx; // the section containing the depth
  }

  for (int i = 0; i < n_depths; i++) { // for each depth...
    double d = depths[i];
    int k = depth_section[i];

    double z0 = elbows[k]; // elbow above
    double frac = (d - z0) / (elbows[k + 1] - z0); // check how far d is into the section

    NumericVector age_vec(n_rows);
    for (int r = 0; r < n_rows; r++) { // now find the age for each MCMC iteration
      double age0 = theta(r, k);
      age_vec[r] = age0 + frac * (theta(r, k + 1) - age0);
    }

    std::sort(age_vec.begin(), age_vec.end()); // place in order

    double qmin = quantileR(age_vec, lower);
    double qmax = quantileR(age_vec, upper);
    double median = quantileR(age_vec, 0.5);
    double mean = std::accumulate(age_vec.begin(), age_vec.end(), 0.0) / n_rows;

    quantiles(i, 0) = depths[i];
    quantiles(i, 1) = qmin;
    quantiles(i, 2) = qmax;
    quantiles(i, 3) = median;
    quantiles(i, 4) = mean;
  }

  colnames(quantiles) = CharacterVector::create("depths", "qmin", "qmax", "median", "mean");
  return quantiles;
}



// [[Rcpp::export]]
NumericMatrix depths_ageranges_hiatus(
  const NumericVector& depths,
  const NumericMatrix& out,
  const NumericVector& elbows, // elbow depths: length K+1
  const NumericVector& hiatus_depths,      // length H
  const NumericMatrix& slopes_above,       // n_rows × H
  const NumericMatrix& slopes_below,       // n_rows × H
  const NumericMatrix& elbow_above_hiatus, // n_rows × H
  const NumericMatrix& elbow_below_hiatus, // n_rows × H
  int n_rows = 4000,
  double prob = 0.95) {

  int n_depths = depths.size(); // number of depths to calculate
  int K = elbows.size() - 1; // number of sections
  double thick = elbows[1] - elbows[0]; // length of sections

  NumericMatrix theta(n_rows, K + 1); // age
  NumericMatrix quantiles(n_depths, 5); // depths, min95%, max95%, median, mean
  double lower = (1.0 - prob) / 2.0; // quantile 0.025 if prob=0.95
  double upper = 1.0 - lower; // quantile 0.975 if prob=0.95

  for (int r = 0; r < n_rows; r++) { // find elbow ages
    theta(r, 0) = out(r, 0); // core-top age for iteration r
    for (int k = 1; k <= K; k++)
      theta(r, k) = theta(r, k - 1) + out(r, k) * thick;
  }

  int H = hiatus_depths.size();
  std::vector<int> hiatus_section(H);

  // hiatus_section[h] = max(which(elbows < hiatus_depths[h]))
  for (int h = 0; h < H; h++) {
    double hd = hiatus_depths[h];
    int idx = std::upper_bound(elbows.begin(), elbows.end(), hd)
      - elbows.begin() - 1;
    if (idx < 0) idx = 0;
    if (idx >= K) idx = K - 1;
    hiatus_section[h] = idx;
  }

  // precompute section index k for each depth
  std::vector<int> depth_section(n_depths);
  for (int i = 0; i < n_depths; i++) {
    double d = depths[i];
    int idx = std::upper_bound(elbows.begin(), elbows.end(), d)
      - elbows.begin() - 1;

    if (idx < 0)   idx = 0;
    if (idx >= K)  idx = K - 1;
    depth_section[i] = idx; // the section containing the depth
  }

  for (int i = 0; i < n_depths; i++) { // for each depth...

    double d = depths[i];
    int k = depth_section[i];
    double z0 = elbows[k]; // elbow above
    double frac = (d - z0) / (elbows[k + 1] - z0); // check how far d is into the section
    NumericVector age_vec(n_rows);

    for (int r = 0; r < n_rows; r++) {
      double age;
      bool in_hiatus = false;

      for (int h = 0; h < H; h++) {
        double z0 = elbows[k];
        double z1 = elbows[k + 1];

        if (k == hiatus_section[h] && d > z0 && d <= z1) {
          double hd = hiatus_depths[h];

          if (d > hd) {
            // BELOW hiatus
            age = elbow_below_hiatus(r,h)
              + slopes_below(r,h) * (d - z1);
          } else {
            // ABOVE hiatus
            age = elbow_above_hiatus(r,h)
              + slopes_above(r,h) * (d - z0);
            }

          in_hiatus = true;
          break;
        }
      }

      // --- normal interpolation if no hiatus applies ---
      if (!in_hiatus) {
        double age0 = theta(r, k);
        age = age0 + frac * (theta(r, k + 1) - age0);
      }

      age_vec[r] = age;
    }

    std::sort(age_vec.begin(), age_vec.end()); // place in order

    double qmin = quantileR(age_vec, lower);
    double qmax = quantileR(age_vec, upper);
    double median = quantileR(age_vec, 0.5);
    double mean = std::accumulate(age_vec.begin(), age_vec.end(), 0.0) / n_rows;

    quantiles(i, 0) = depths[i];
    quantiles(i, 1) = qmin;
    quantiles(i, 2) = qmax;
    quantiles(i, 3) = median;
    quantiles(i, 4) = mean;
  }

  colnames(quantiles) = CharacterVector::create("depths", "qmin", "qmax", "median", "mean");
  return quantiles;
}



// [[Rcpp::export]]
List depths_agegrid(
  const NumericVector& depths,
  const NumericMatrix& out,
  const NumericVector& elbows, // elbow depths: length K+1
  int hist_n = 512,
  double min_age = 0.0,
  double max_age = 55000.0,
  int n_rows = 4000,
  double prob = 0.95) {

  int n_depths = depths.size(); // number of depths to calculate
  int K = elbows.size() - 1; // number of sections
  double thick = elbows[1] - elbows[0]; // length of sections

  NumericMatrix theta(n_rows, K + 1); // age
  int n_bins = hist_n - 1;
  NumericVector breaks(hist_n);
  NumericMatrix densities(n_depths, n_bins);
  double bin_width = (max_age - min_age) / n_bins;
  for (int b = 0; b < hist_n; b++) {
    breaks[b] = min_age + b * bin_width;
  }

  for (int r = 0; r < n_rows; r++) { // find elbow ages
    theta(r, 0) = out(r, 0); // core-top age for iteration r
    for (int k = 1; k <= K; k++)
      theta(r, k) = theta(r, k - 1) + out(r, k) * thick;
  }

  // precompute section index k for each depth
  std::vector<int> depth_section(n_depths);
  for (int i = 0; i < n_depths; i++) {
    double d = depths[i];
    int idx = std::upper_bound(elbows.begin(), elbows.end(), d)
      - elbows.begin() - 1;

    if (idx < 0)   idx = 0;
    if (idx >= K)  idx = K - 1;
    depth_section[i] = idx; // the section containing the depth
  }

  for (int i = 0; i < n_depths; i++) { // for each depth...
    double d = depths[i];
    int k = depth_section[i];

    double z0 = elbows[k]; // elbow above
    double frac = (d - z0) / (elbows[k + 1] - z0); // check how far d is into the section

    NumericVector age_vec(n_rows);
    NumericVector dens(n_bins);

    for (int r = 0; r < n_rows; r++) { // now find the age for each MCMC iteration
      double age0 = theta(r, k);
      age_vec[r] = age0 + frac * (theta(r, k + 1) - age0);
    }

    for (int r = 0; r < n_rows; r++) {
      int bin = floor((age_vec[r] - min_age) / bin_width);
      if (bin < 0) continue;
      if (bin >= n_bins) bin = n_bins -1;
      dens[bin] += 1.0;
    }

    densities(i, _) = dens / n_rows;
  }

  return List::create(
    Named("density")   = densities,
    Named("breaks")    = breaks,
    Named("depths")    = depths);
}



// [[Rcpp::export]]
List depths_agegrid_hiatus(
  const NumericVector& depths,
  const NumericMatrix& out,
  const NumericVector& elbows, // elbow depths: length K+1
  const NumericVector& hiatus_depths,      // length H
  const NumericMatrix& slopes_above,       // n_rows × H
  const NumericMatrix& slopes_below,       // n_rows × H
  const NumericMatrix& elbow_above_hiatus, // n_rows × H
  const NumericMatrix& elbow_below_hiatus, // n_rows × H
  int hist_n = 512,
  double min_age = 0.0,
  double max_age = 55000.0,
  int n_rows = 4000,
  double prob = 0.95) {

  int n_depths = depths.size(); // number of depths to calculate
  int K = elbows.size() - 1; // number of sections
  double thick = elbows[1] - elbows[0]; // length of sections

  NumericMatrix theta(n_rows, K + 1); // age
  int n_bins = hist_n - 1;
  NumericVector breaks(hist_n);
  NumericMatrix densities(n_depths, n_bins);
  double bin_width = (max_age - min_age) / n_bins;
  for (int b = 0; b < hist_n; b++) {
    breaks[b] = min_age + b * bin_width;
  }

  for (int r = 0; r < n_rows; r++) { // find elbow ages
    theta(r, 0) = out(r, 0); // core-top age for iteration r
    for (int k = 1; k <= K; k++)
      theta(r, k) = theta(r, k - 1) + out(r, k) * thick;
  }

  int H = hiatus_depths.size();
  std::vector<int> hiatus_section(H);

  // hiatus_section[h] = max(which(elbows < hiatus_depths[h]))
  for (int h = 0; h < H; h++) {
    double hd = hiatus_depths[h];
    int idx = std::upper_bound(elbows.begin(), elbows.end(), hd)
      - elbows.begin() - 1;
    if (idx < 0) idx = 0;
    if (idx >= K) idx = K - 1;
    hiatus_section[h] = idx;
  }

  // precompute section index k for each depth
  std::vector<int> depth_section(n_depths);
  for (int i = 0; i < n_depths; i++) {
    double d = depths[i];
    int idx = std::upper_bound(elbows.begin(), elbows.end(), d)
      - elbows.begin() - 1;

    if (idx < 0)   idx = 0;
    if (idx >= K)  idx = K - 1;
    depth_section[i] = idx; // the section containing the depth
  }

  for (int i = 0; i < n_depths; i++) { // for each depth...
    double d = depths[i];
    int k = depth_section[i];
    double z0 = elbows[k]; // elbow above
    double frac = (d - z0) / (elbows[k + 1] - z0); // check how far d is into the section
    NumericVector age_vec(n_rows);
    NumericVector dens(n_bins);

    for (int r = 0; r < n_rows; r++) {
      double age;
      bool in_hiatus = false;

      for (int h = 0; h < H; h++) {
        double z0 = elbows[k];
        double z1 = elbows[k + 1];

        if (k == hiatus_section[h] && d > z0 && d <= z1) {
          double hd = hiatus_depths[h];

          if (d > hd) {
            // BELOW hiatus
            age = elbow_below_hiatus(r,h)
              + slopes_below(r,h) * (d - z1);
          } else {
            // ABOVE hiatus
            age = elbow_above_hiatus(r,h)
              + slopes_above(r,h) * (d - z0);
            }

          in_hiatus = true;
          break;
        }
      }

      // --- normal interpolation if no hiatus applies ---
      if (!in_hiatus) {
        double age0 = theta(r, k);
        age = age0 + frac * (theta(r, k + 1) - age0);
      }

      age_vec[r] = age;
    }

    for (int r = 0; r < n_rows; r++) {
      int bin = floor((age_vec[r] - min_age) / bin_width);
      if (bin < 0) continue;
      if (bin >= n_bins) bin = n_bins - 1;
      dens[bin] += 1.0;
    }

    densities(i, _) = dens / n_rows;
  }

  return List::create(
      Named("density") = densities,
      Named("breaks") = breaks,
      Named("depths") = depths);
}

