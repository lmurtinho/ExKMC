#ifndef BEST_CUT
#define BEST_CUT

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

struct Cut {
  int dim;
  double cut;
  double cost;
};

void best_cut_dim(double *data, double *centers, double *distances,
                  int *dist_order, double *cuts, double *costs,
                  int *data_order, int *centers_order, int n, int d,
                  int k, int n_cuts, int dim, double *ans) {
  int i, j, c, cur_c, new_c, ix;
  int idx_cut = 0;
  int idx_data = 0;
  int idx_centers = 0;
  int data_in_left;
  int center_in_left;
  double nxt_cut;
  double best_cut;
  double cur_cost = 0;
  double best_cost;
  double old_data_cost;
  // double sums;

  printf("I'm in!\n");
  //
  // sums = 0;
  // for (i = 0; i < n_cuts; i++) {
  //   sums += cuts[i];
  // }
  // printf("Sum of %d possible cuts for dimension %d: %.4f\n", n_cuts, dim, sums);
  //
  // for (i = 0; i < n_cuts; i++) {
  //   printf("%.4f\t", cuts[i]);
  // }
  // printf("\n");

  // current costs are the best costs (no cuts yet)
  double *cur_costs = (double *)malloc(sizeof(double) * n);
  for (i = 0; i < n; i++) {
    cur_costs[i] = costs[i];
    cur_cost += costs[i];
  }
  // printf("Initial cost: %.4f\n", cur_cost);

  // all data starts to the right
  int *left_data_mask = (int *)malloc(sizeof(int) * n);
  for (i = 0; i < n; i++) {
    left_data_mask[i] = 0;
  }
  int *left_centers_mask = (int *)malloc(sizeof(int) * k);
  for (i = 0; i < k; i++) {
    left_centers_mask[i] = 0;
  }

  // printf("masks done\n");

  // CHECK: current centers are first in each row of centers_order
  int *cur_centers = (int *)malloc(sizeof(int) * n);
  for (i = 0; i < n; i++) {
    cur_centers[i] = dist_order[i * k];
  }

  // for (i = 0; i < n; i++) {
  //   if ( (cur_centers[i] < 0) || (cur_centers[i] >= k) ) {
  //     printf("Incorrect center value for %d: %d\n", i, cur_centers[i]);
  //   }
  // }

  // printf("Got current centers!\n");

  // the first cut moves the first center to the left
  nxt_cut = centers[centers_order[0]*d + dim];
  while (cuts[idx_cut] != nxt_cut) {
    idx_cut++;
  }
  // idx_cut++;

  // printf("Got first cut!\n");

  // for (i = 0; i < k; i++) {
  //   printf("%d\t", left_centers_mask[i]);
  // }
  // printf("\n");

  // set indices of data and center as first index in which the element is to
  // the right of the cut
  while (data[data_order[idx_data]*d + dim] <= nxt_cut) {
    idx_data++;
  }
  while (centers[centers_order[idx_centers]*d + dim] <= nxt_cut) {
    idx_centers++;
  }

  // printf("Data and center indices are set!\n");

  // redefine masks
  for (i=0; i < idx_data; i++) {
    left_data_mask[data_order[i]] = 1;
  }
  for (i = 0; i < idx_centers; i++) {
    left_centers_mask[centers_order[i]] = 1;
  }

  // for (i = 0; i < k; i++) {
  //   printf("%d\t", left_centers_mask[i]);
  // }
  // printf("\n");
  //
  // printf("Masks are defined!\n");

  // reassign data points separated from their centers
  for (i = 0; i < n; i++) {
    // printf("Dealing with point %d \n", i);
    cur_c = cur_centers[i];
    if (left_data_mask[i] != left_centers_mask[cur_c]) {
      // printf("Need to reassign point %d \n", i);
      old_data_cost = cur_costs[i];
      // printf("Old cost = %.4f \n", old_data_cost);
      data_in_left = left_data_mask[i];
      // printf("Is data to the left? %d\n", data_in_left);
      j = 0;
      while (data_in_left != left_centers_mask[centers_order[j]]) {
        // printf("Is center %d to the left? %d\n", centers_order[j], left_centers_mask[centers_order[j]]);
        j++;
      }
      c = centers_order[j];
      // printf("Is center %d to the left? %d\n", c, left_centers_mask[c]);
      cur_centers[i] = c;
      cur_costs[i] = distances[i*k + c];
      // printf("New cost for datum: %.4f\n", cur_costs[i]);
      cur_cost += (cur_costs[i] - old_data_cost);
      // printf("New cost for data: %.4f\n", cur_cost);
    }
  }

  // printf("Initialization done!\n");

  // store initial cut and cost as best ones
  best_cut = nxt_cut;
  best_cost = cur_cost;

  // printf("cut: %.4f\tcost: %.4f\n", nxt_cut, cur_cost);

  while ( (idx_cut < n_cuts) && (idx_centers < k) ) {
    // printf("cut %d of %d\n", idx_cut, n_cuts);
    // printf("center %d of %d\n", idx_centers, k);
    // printf("datum %d of %d\n", idx_data, n);
    nxt_cut = cuts[idx_cut];
    // move data points to the left and assign them to best center to the left
    while ( (idx_data < n) && (data[data_order[idx_data] * d + dim] <= nxt_cut) ) {
      i = data_order[idx_data];
      old_data_cost = cur_costs[i];
      left_data_mask[i] = 1;
      j = 0;
      while (!left_centers_mask[dist_order[i*k + j]]) {
        j++;
      }
      c = dist_order[i*k + j];
      cur_centers[i] = c;
      cur_costs[i] = distances[i*k + c];
      cur_cost += (cur_costs[i] - old_data_cost);
      idx_data++;
    }

    // move centers to the left and reassign data points
    while ( (idx_centers < k) && (centers[centers_order[idx_centers]*d + dim] <= nxt_cut) ) {
      c = centers_order[idx_centers];
      left_centers_mask[c] = 1;
      for (i = 0; i < n; i++) {
        // if datum is to the left, assign it to center that moved to the left
        // if preferable to current center
        if (left_data_mask[i]) {
          old_data_cost = cur_costs[i];
          if (distances[i*k + c] < old_data_cost) {
            cur_costs[i] = distances[i*k + c];
            cur_centers[i] = c;
            cur_cost += (cur_costs[i] - old_data_cost);
          }
        }
        // if datum is to the right and its current center moved left, find
        // new best center
        else if (cur_centers[i] == c) {
          old_data_cost = cur_costs[i];
          j = 0;
          while (left_centers_mask[dist_order[i*k + j]]) {
            j++;
          }
          new_c = dist_order[i*k + j];
          cur_costs[i] = distances[i*k + new_c];
          cur_cost += (cur_costs[i] - old_data_cost);
        }
      }
      idx_centers++;
    }
    // replace best cut if current cost is smaller than best cost
    // printf("cut: %.4e\tcost: %.4e\n", nxt_cut, cur_cost);
    // printf("costs: %.4f\t%.4f\n", cur_cost, best_cost);
    if (cur_cost < best_cost) {
      best_cut = nxt_cut;
      best_cost = cur_cost;
      // ans.cut = best_cut;
      // ans.cost = best_cost;
      // printf("ans: %.4f\t%.4f\n", ans.cut, ans.cost);
      // printf("best cut: %.4e\tbest cost: %.4e\n", best_cut, best_cost);
    }
    idx_cut++;
  }

  // printf("Done!\n");

  ans[0] = best_cut;
  ans[1] = best_cost;

  free(cur_costs);
  free(left_data_mask);
  free(left_centers_mask);
  free(cur_centers);
  // printf("ready to return!\n");

  return;
};

void best_cut_single_dim(double *data, double *centers, double *distances,
                         int *dist_order, double *cuts, double *costs,
                         int *data_order, int *centers_order, int n, int k,
                         int n_cuts, double *ans) {
  // printf("start: %d, %d, %d", n, k, n_cuts);
  int i, j, c, cur_c, new_c, ix;
  int idx_cut = 0;
  int idx_data = 0;
  int idx_centers = 0;
  int data_in_left;
  int center_in_left;
  double nxt_cut;
  double best_cut;
  double cur_cost = 0;
  double best_cost;
  double old_data_cost;
  // double sums;

  // printf("I'm in!\n");
  //
  // sums = 0;
  // for (i = 0; i < n_cuts; i++) {
  //   sums += cuts[i];
  // }
  // printf("Sum of %d possible cuts for dimension %d: %.4f\n", n_cuts, dim, sums);
  //
  // for (i = 0; i < n_cuts; i++) {
  //   printf("%.4f\t", cuts[i]);
  // }
  // printf("\n");

  // current costs are the best costs (no cuts yet)
  double *cur_costs = (double *)malloc(sizeof(double) * n);
  for (i = 0; i < n; i++) {
    cur_costs[i] = costs[i];
    cur_cost += costs[i];
  }
  // printf("Initial cost: %.4f\n", cur_cost);

  // all data starts to the right
  int *left_data_mask = (int *)malloc(sizeof(int) * n);
  for (i = 0; i < n; i++) {
    left_data_mask[i] = 0;
  }
  int *left_centers_mask = (int *)malloc(sizeof(int) * k);
  for (i = 0; i < k; i++) {
    left_centers_mask[i] = 0;
  }

  // printf("masks done\n");

  // CHECK: current centers are first in each row of centers_order
  int *cur_centers = (int *)malloc(sizeof(int) * n);
  for (i = 0; i < n; i++) {
    cur_centers[i] = dist_order[i * k];
  }

  // for (i = 0; i < n; i++) {
  //   if ( (cur_centers[i] < 0) || (cur_centers[i] >= k) ) {
  //     printf("Incorrect center value for %d: %d\n", i, cur_centers[i]);
  //   }
  // }

  // printf("Got current centers!\n");

  // the first cut moves the first center to the left
  nxt_cut = centers[centers_order[0]];
  while (cuts[idx_cut] != nxt_cut) {
    idx_cut++;
  }
  // idx_cut++;

  // printf("Got first cut!\n");

  // for (i = 0; i < k; i++) {
  //   printf("%d\t", left_centers_mask[i]);
  // }
  // printf("\n");

  // set indices of data and center as first index in which the element is to
  // the right of the cut
  while (data[data_order[idx_data]] <= nxt_cut) {
    idx_data++;
  }
  while (centers[centers_order[idx_centers]] <= nxt_cut) {
    idx_centers++;
  }

  // printf("Data and center indices are set!\n");

  // redefine masks
  for (i=0; i < idx_data; i++) {
    left_data_mask[data_order[i]] = 1;
  }
  for (i = 0; i < idx_centers; i++) {
    left_centers_mask[centers_order[i]] = 1;
  }

  // for (i = 0; i < k; i++) {
  //   printf("%d\t", left_centers_mask[i]);
  // }
  // printf("\n");
  //
  // printf("Masks are defined!\n");

  // reassign data points separated from their centers
  for (i = 0; i < n; i++) {
    // printf("Dealing with point %d \n", i);
    cur_c = cur_centers[i];
    if (left_data_mask[i] != left_centers_mask[cur_c]) {
      // printf("Need to reassign point %d \n", i);
      old_data_cost = cur_costs[i];
      // printf("Old cost = %.4f \n", old_data_cost);
      data_in_left = left_data_mask[i];
      // printf("Is data to the left? %d\n", data_in_left);
      j = 0;
      while (data_in_left != left_centers_mask[centers_order[j]]) {
        // printf("Is center %d to the left? %d\n", centers_order[j], left_centers_mask[centers_order[j]]);
        j++;
      }
      c = centers_order[j];
      // printf("Is center %d to the left? %d\n", c, left_centers_mask[c]);
      cur_centers[i] = c;
      cur_costs[i] = distances[i*k + c];
      // printf("New cost for datum: %.4f\n", cur_costs[i]);
      cur_cost += (cur_costs[i] - old_data_cost);
      // printf("New cost for data: %.4f\n", cur_cost);
    }
  }

  // printf("Initialization done!\n");

  // store initial cut and cost as best ones
  best_cut = nxt_cut;
  best_cost = cur_cost;

  // printf("cut: %.4f\tcost: %.4f\n", nxt_cut, cur_cost);

  while ( (idx_cut < n_cuts) && (idx_centers < k) ) {
    // printf("cut %d of %d\n", idx_cut, n_cuts);
    // printf("center %d of %d\n", idx_centers, k);
    // printf("datum %d of %d\n", idx_data, n);
    nxt_cut = cuts[idx_cut];
    // move data points to the left and assign them to best center to the left
    while ( (idx_data < n) && (data[data_order[idx_data]] <= nxt_cut) ) {
      i = data_order[idx_data];
      old_data_cost = cur_costs[i];
      left_data_mask[i] = 1;
      j = 0;
      while (!left_centers_mask[dist_order[i*k + j]]) {
        j++;
      }
      c = dist_order[i*k + j];
      cur_centers[i] = c;
      cur_costs[i] = distances[i*k + c];
      cur_cost += (cur_costs[i] - old_data_cost);
      idx_data++;
    }

    // move centers to the left and reassign data points
    while ( (idx_centers < k) && (centers[centers_order[idx_centers]] <= nxt_cut) ) {
      c = centers_order[idx_centers];
      left_centers_mask[c] = 1;
      for (i = 0; i < n; i++) {
        // if datum is to the left, assign it to center that moved to the left
        // if preferable to current center
        if (left_data_mask[i]) {
          old_data_cost = cur_costs[i];
          if (distances[i*k + c] < old_data_cost) {
            cur_costs[i] = distances[i*k + c];
            cur_centers[i] = c;
            cur_cost += (cur_costs[i] - old_data_cost);
          }
        }
        // if datum is to the right and its current center moved left, find
        // new best center
        else if (cur_centers[i] == c) {
          old_data_cost = cur_costs[i];
          j = 0;
          while (left_centers_mask[dist_order[i*k + j]]) {
            j++;
          }
          new_c = dist_order[i*k + j];
          cur_costs[i] = distances[i*k + new_c];
          cur_cost += (cur_costs[i] - old_data_cost);
          cur_centers[i] = new_c;
        }
      }
      idx_centers++;
    }
    // replace best cut if current cost is smaller than best cost
    // printf("cut: %.4e\tcost: %.4e\n", nxt_cut, cur_cost);
    // printf("costs: %.4f\t%.4f\n", cur_cost, best_cost);
    if (cur_cost < best_cost) {
      best_cut = nxt_cut;
      best_cost = cur_cost;
      // ans.cut = best_cut;
      // ans.cost = best_cost;
      // printf("ans: %.4f\t%.4f\n", ans.cut, ans.cost);
      // printf("best cut: %.4e\tbest cost: %.4e\n", best_cut, best_cost);
    }
    idx_cut++;
  }

  // printf("Done!\n");

  ans[0] = best_cut;
  ans[1] = best_cost;

  free(cur_costs);
  free(left_data_mask);
  free(left_centers_mask);
  free(cur_centers);
  printf("\tready to return!\n");

  return;
};

void best_cut_single_dim_opt(double *data, double *centers, double *distances,
                             int *dist_order, int *data_order,
                             int *centers_order, int n, int k, double *ans) {

  int i, j, c, cur_c, new_c, ix, ic;
  int idx_data = 0;
  int idx_centers = 0;
  int data_in_left;
  int center_in_left;
  double nxt_cut;
  double best_cut;
  double cur_cost = 0;
  double best_cost;
  double old_data_cost;
  double max_cut;

  // max cut is largest value of center
  ic = centers_order[k-1];
  max_cut = centers[ic];

  // current costs are the best costs (no cuts yet)
  double *cur_costs = (double *)malloc(sizeof(double) * n);
  int *cur_centers = (int *)malloc(sizeof(int) * n);
  for (i = 0; i < n; i++) {
    c = dist_order[i * k];
    cur_centers[i] = c;
    cur_costs[i] = distances[i*k + c];
    cur_cost += cur_costs[i];
  }

  // all data starts to the right
  int *left_data_mask = (int *)malloc(sizeof(int) * n);
  for (i = 0; i < n; i++) {
    left_data_mask[i] = 0;
  }
  int *left_centers_mask = (int *)malloc(sizeof(int) * k);
  for (i = 0; i < k; i++) {
    left_centers_mask[i] = 0;
  }


  // the first cut moves the first center to the left
  c = centers_order[0];
  nxt_cut = centers[c];

  // set indices of data and center as first index in which the element is to
  // the right of the cut
  ix = data_order[idx_data];
  while ( (data[ix] <= nxt_cut) && idx_data < n) {
    idx_data++;
    if (idx_data < n) {
      ix = data_order[idx_data];
    }
  }

  while ( (centers[c] <= nxt_cut) && (idx_centers < k) ) {
    idx_centers++;
    if (idx_centers < k) {
      c = centers_order[idx_centers];
      }
    }

  // redefine masks
  for (i=0; i < idx_data; i++) {
    ix = data_order[i];
    left_data_mask[ix] = 1;
  }
  for (i = 0; i < idx_centers; i++) {
    ic = centers_order[i];
    left_centers_mask[ic] = 1;
  }

  // if all centers moved to the left after the
  // first cut, return (no cut is possible)
  if (idx_centers == k) {
    ans[0] = -1;
    ans[1] = INFINITY;
    free(cur_costs);
    free(left_data_mask);
    free(left_centers_mask);
    free(cur_centers);
    return;
  }

  // reassign data points separated from their centers
  for (i = 0; i < n; i++) {
    c = cur_centers[i];
    if (left_data_mask[i] != left_centers_mask[c]) {
      old_data_cost = cur_costs[i];
      data_in_left = left_data_mask[i];
      j = 0;
      while (data_in_left != left_centers_mask[c]) {
        j++;
        c = dist_order[i*k + j];
      }
      cur_centers[i] = c;
      cur_costs[i] = distances[i*k + c];
      cur_cost += (cur_costs[i] - old_data_cost);
    }
  }

  // store initial cut and cost as best ones
  best_cut = nxt_cut;
  best_cost = cur_cost;

  while ( (idx_data < n) && (idx_centers < k) ) {
    // find next cut and check if feasible
    ix = data_order[idx_data];
    ic = centers_order[idx_centers];
    if (data[ix] < centers[ic]) {
      nxt_cut = data[ix];
    }
    else {
      nxt_cut = centers[ic];
    }
    if (nxt_cut >= max_cut) {
      break;
    }

    // move data points to the left and assign them to best center to the left
    while ( (idx_data < n) && (data[ix] <= nxt_cut) ) {
      old_data_cost = cur_costs[ix];
      left_data_mask[ix] = 1;
      j = 0;
      c = dist_order[ix*k];
      while (!left_centers_mask[c]) {
        j++;
        c = dist_order[ix*k + j];
      }
      cur_centers[ix] = c;
      cur_costs[ix] = distances[ix*k + c];
      cur_cost += (cur_costs[ix] - old_data_cost);
      idx_data++;
      if (idx_data < n) {
          ix = data_order[idx_data];
      }
    }

    // move centers to the left and reassign data points
    while ( (idx_centers < k) && (centers[ic] <= nxt_cut) ) {
      left_centers_mask[ic] = 1;
      for (i = 0; i < n; i++) {
        old_data_cost = cur_costs[i];
        // if datum is to the left, assign it to center that moved to the left
        // if it's best in left
        if (left_data_mask[i]) {
          if (distances[i*k + ic] < old_data_cost) {
            cur_centers[i] = ic;
            cur_costs[i] = distances[i*k + ic];
            cur_cost += (cur_costs[i] - old_data_cost);
          }
        }
        // if datum is to the right and its current center moved left, find
        // new best center
        else if (cur_centers[i] == ic) {
          j = 0;
          c = dist_order[i*k];
          while (left_centers_mask[c]) {
            j++;
            c = dist_order[i*k + j];
          }
          cur_costs[i] = distances[i*k + c];
          cur_cost += (cur_costs[i] - old_data_cost);
          cur_centers[i] = c;
        }
      }
      idx_centers++;
      if (idx_centers < k) {
        ic = centers_order[idx_centers];
      }
    }
    // replace best cut if current cost is smaller than best cost
    if (cur_cost < best_cost) {
      best_cut = nxt_cut;
      best_cost = cur_cost;
    }
  }

  // return best cut and cost

  ans[0] = best_cut;
  ans[1] = best_cost;

  free(cur_costs);
  free(left_data_mask);
  free(left_centers_mask);
  free(cur_centers);

  return;
};

void best_cut_single_dim_opt2(double *data, double *centers, double *distances,
                             int *dist_order, int *data_order,
                             int *centers_order, int n, int k, double *ans) {

  int i, j, c, cur_c, new_c, ix, ic;
  int idx_data = 0;
  int idx_centers = 0;
  int data_in_left;
  int center_in_left;
  double nxt_cut;
  double best_cut;
  double cur_cost = 0;
  double best_cost;
  double old_data_cost;
  double max_cut;

  // max cut is largest value of center
  ic = centers_order[k-1];
  max_cut = centers[ic];

  // current costs are the best costs (no cuts yet)
  double *cur_costs = (double *)malloc(sizeof(double) * n);
  int *cur_centers = (int *)malloc(sizeof(int) * n);
  for (i = 0; i < n; i++) {
    c = dist_order[i * k];
    cur_centers[i] = c;
    cur_costs[i] = distances[i*k + c];
    cur_cost += cur_costs[i];
  }

  // all data starts to the right
  int *left_data_mask = (int *)malloc(sizeof(int) * n);
  for (i = 0; i < n; i++) {
    left_data_mask[i] = 0;
  }
  int *left_centers_mask = (int *)malloc(sizeof(int) * k);
  for (i = 0; i < k; i++) {
    left_centers_mask[i] = 0;
  }


  // the first cut moves the first center to the left
  c = centers_order[0];
  nxt_cut = centers[c];

  // CHANGED: added vector to keep track of best center to the left
  // initially, the best center to the left is the only one there
  int *best_in_left = (int *)malloc(sizeof(int) * n);
  for (i = 0; i < n; i++) {
    best_in_left[i] = c;
  }

  // set indices of data and center as first index in which the element is to
  // the right of the cut
  ix = data_order[idx_data];
  while ( (data[ix] <= nxt_cut) && idx_data < n) {
    idx_data++;
    if (idx_data < n) {
      ix = data_order[idx_data];
    }
  }

  while ( (centers[c] <= nxt_cut) && (idx_centers < k) ) {
    idx_centers++;
    if (idx_centers < k) {
      c = centers_order[idx_centers];
      // CHANGED: if a center is moved to the left duting initialization,
      // check if it's the best center to the left for each datum
      for (i = 0; i < n; i++) {
        cur_c = best_in_left[i];
        if (distances[i*k + cur_c] > distances[i*k + c]) {
          best_in_left[i] = c;
        }
      }
    }
  }

  // redefine masks
  for (i=0; i < idx_data; i++) {
    ix = data_order[i];
    left_data_mask[ix] = 1;
  }
  for (i = 0; i < idx_centers; i++) {
    ic = centers_order[i];
    left_centers_mask[ic] = 1;
  }

  // if all centers moved to the left after the
  // first cut, return (no cut is possible)
  if (idx_centers == k) {
    ans[0] = -1;
    ans[1] = INFINITY;
    free(cur_costs);
    free(left_data_mask);
    free(left_centers_mask);
    free(cur_centers);
    return;
  }

  // reassign data points separated from their centers
  // for (i = 0; i < n; i++) {
  //   c = cur_centers[i];
  //   if (left_data_mask[i] != left_centers_mask[c]) {
  //     old_data_cost = cur_costs[i];
  //     data_in_left = left_data_mask[i];
  //     j = 0;
  //     while (data_in_left != left_centers_mask[c]) {
  //       j++;
  //       c = dist_order[i*k + j];
  //     }
  //     cur_centers[i] = c;
  //     cur_costs[i] = distances[i*k + c];
  //     cur_cost += (cur_costs[i] - old_data_cost);
  //   }
  // }

  // reassign data separated from their centers
  for (i = 0; i < n; i++) {
    cur_c = cur_centers[i];
    // if datum is to the left, best center is in best_in_left vector
    if (left_data_mask[i]) {
      c = best_in_left[i];
    }
    // if datum is to the right, find center to the right closest to datum
    else {
      j = 0;
      c = dist_order[i*k + j];
      while (left_centers_mask[c]) {
        j++;
        c = dist_order[i*k + j];
      }
    }
    // if best center is different than current center, update cost
    if (c != cur_c) {
      old_data_cost = cur_costs[i];
      cur_costs[i] = distances[i*k + c];
      cur_cost += (cur_costs[i] - old_data_cost);
      cur_centers[i] = c;
    }
  }

  // store initial cut and cost as best ones
  best_cut = nxt_cut;
  best_cost = cur_cost;

  while ( (idx_data < n) && (idx_centers < k) ) {
    // find next cut and check if feasible
    ix = data_order[idx_data];
    ic = centers_order[idx_centers];
    if (data[ix] < centers[ic]) {
      nxt_cut = data[ix];
    }
    else {
      nxt_cut = centers[ic];
    }
    if (nxt_cut >= max_cut) {
      break;
    }

    // move data points to the left and assign them to best center to the left
    while ( (idx_data < n) && (data[ix] <= nxt_cut) ) {
      old_data_cost = cur_costs[ix];
      left_data_mask[ix] = 1;
      // CHANGED: find best center to the left via best_in_left vector
      c = best_in_left[ix];
      cur_centers[ix] = c;
      cur_costs[ix] = distances[ix*k + c];
      cur_cost += (cur_costs[ix] - old_data_cost);
      idx_data++;
      if (idx_data < n) {
          ix = data_order[idx_data];
      }
    }

    // move centers to the left and reassign data points
    while ( (idx_centers < k) && (centers[ic] <= nxt_cut) ) {
      left_centers_mask[ic] = 1;
      for (i = 0; i < n; i++) {
        old_data_cost = cur_costs[i];
        // CHANGED: update best_in_left vector
        cur_c = best_in_left[i];
        if (distances[i*k + ic] < distances[i*k + cur_c]) {
          best_in_left[i] = ic;
        }
        // if datum is to the left, assign it to center that moved to the left
        // if it's best in left
        if (left_data_mask[i]) {
          if (best_in_left[i] == ic) {
            cur_centers[i] = ic;
            cur_costs[i] = distances[i*k + ic];
            cur_cost += (cur_costs[i] - old_data_cost);
          }
        }
        // if datum is to the right and its current center moved left, find
        // new best center
        else if (cur_centers[i] == ic) {
          j = 0;
          c = dist_order[i*k];
          while (left_centers_mask[c]) {
            j++;
            c = dist_order[i*k + j];
          }
          cur_costs[i] = distances[i*k + c];
          cur_cost += (cur_costs[i] - old_data_cost);
          cur_centers[i] = c;
        }
      }
      idx_centers++;
      if (idx_centers < k) {
        ic = centers_order[idx_centers];
      }
    }
    // replace best cut if current cost is smaller than best cost
    if (cur_cost < best_cost) {
      best_cut = nxt_cut;
      best_cost = cur_cost;
    }
  }

  // return best cut and cost

  ans[0] = best_cut;
  ans[1] = best_cost;

  free(cur_costs);
  free(left_data_mask);
  free(left_centers_mask);
  free(cur_centers);

  return;
};

double norm2squared(double *u, double *v, int d) {
  int i;
  double ans = 0;
  double val;
  for (i = 0; i < d; i++) {
    val = u[i] - v[i];
    val *= val;
    ans += val;
  }
  return ans;
}

double *get_distances(double *data, double *centers, int n, int k, int d) {
  int i, j;
  double *ans = (double *)malloc(sizeof(double) * n * k);
  for (i = 0; i < n; i++) {
    for (j = 0; j < k; j++) {
      ans[i*d + j] = norm2squared(&data[i*d], &centers[j*d], d);
    }
  }
  return ans;
}



#endif
