//
// Demonstration program for hierarchical clustering
// with fastcluster by Daniel Muellner
//
// Line segments are clustered in two directions
//
// Author: Christoph Dalitz, 2018
//

#include <math.h>
#include <string.h>
#include <unistd.h>
#include <cstdint>

#include <string>
#include <vector>
#include <stdio.h>
#include <algorithm>
#include <fstream>
#include <iostream>

#include "fastcluster.h"


// 2D point or vector
class Point {
public:
  double x;
  double y;
  Point(double xx=0.0, double yy=0.0) { x=xx; y=yy; }
  Point(const Point& p) { x=p.x; y=p.y; }
  double norm() { return(sqrt(x*x + y*y)); }
};

// line segment
class Segment {
public:
  Point p1;
  Point p2;
  Point dir;
  Segment(const Point p, const Point q) { p1=p, p2=q; dir=direction(); }
  Segment(const Segment& s) { p1=s.p1, p2=s.p2; dir=direction(); }
private:
  Point direction() {
    Point d(p2.x-p1.x, p2.y-p1.y);
    double n = d.norm();
    if (n > 0.0) {
      d.x /= n; d.y /= n;
    }
    return d;
  }
};

// // line segment distance (cosine dissimilarity)
// double distance(const Segment& s1, const Segment& s2) {
//   double sprod = s1.dir.x*s2.dir.x + s1.dir.y*s2.dir.y;
//   double d = 1 - sprod*sprod;
//   if (d < 0.0)
//     return 0;
//   else
//     return d;
// }

double distance(const double& s1, const double& s2) {
  return fabs(s1-s2);
}

double distance_event(const double& s11, const double& s12, const double& s21, const double& s22) {
  if (s12 < s21)
    return s21 - s12;
  return s11 - s22;
}

int changeMerge(int i, int N){
  if (i < 0) {
    return (i+1) * -1;
  } else {
    return (i - 1) + N;
  }
}

// trying the search query
// search logic
// if bin size is less than cluster length, then go down
// else return the start end point of the current cluster
// dfs on the start and end time query
void dfs(int *merge, int npoints, int start_t, int end_t, int bin_size, int c_node, double* start_events, double* end_events) {
  int root_node = c_node + npoints;
  int start_time = start_events[root_node];
  int end_time = end_events[root_node];
  if(start_time >= end_t || end_time <= start_t) {
    return;
  }
  if(bin_size >= (end_time - start_time + 1)) {
    std::cout << "Cluster: " << root_node << ", Start: " << start_time << ", End: " << end_time << std::endl;
    return;
  }
  if(root_node < npoints) {
    if(start_time < start_t) {
      start_time = start_t;
    }
    if(end_time > end_t) {
      end_time = end_t;
    }
    std::cout << "Cluster-Leaf: " << root_node << ", Start: " << start_time << ", End: " << end_time << std::endl;
    return;
  }
  // this is a compound node
  int left_node_index = changeMerge(merge[c_node], npoints);
  int right_node_index = changeMerge(merge[c_node+npoints-1], npoints);
  if(end_events[left_node_index] > start_events[right_node_index]) {
    std::swap(left_node_index, right_node_index);
  }
  dfs(merge, npoints, start_t, end_t, bin_size, left_node_index - npoints, start_events, end_events);
  dfs(merge, npoints, start_t, end_t, bin_size, right_node_index - npoints, start_events, end_events);
}

int main(int argc, char** argv)
{
  std::string opt_infile;
  int i,j,k,npoints;
  int opt_method = HCLUST_METHOD_SINGLE;

  const char* usagemsg = "Usage: hclust-demo <infile> [-m (single|complete|average|median)]\n";
  for (i=1; i<argc; i++) {
    if (0 == strcmp(argv[i], "-m")) {
      i++;
      if (i<argc) {
        if (0 == strcmp(argv[i], "single"))
            opt_method = HCLUST_METHOD_SINGLE;
        else if (0 == strcmp(argv[i], "complete"))
            opt_method = HCLUST_METHOD_COMPLETE;
        else if (0 == strcmp(argv[i], "average"))
            opt_method = HCLUST_METHOD_AVERAGE;
        else if (0 == strcmp(argv[i], "median"))
            opt_method = HCLUST_METHOD_MEDIAN;
        else {
          fputs(usagemsg, stderr);
          return 1;
        }
      } else {
        fputs(usagemsg, stderr);
        return 1;
      }
    }
    else if (argv[i][0] == '-') {
      fputs(usagemsg, stderr);
      return 1;
    }
    else {
      opt_infile = argv[i];
    }
  }

  // std::vector<double> v = {2, 8, 0, 4, 1, 9, 9, 0, 0, 10};
  std::vector<double> v;// = {1,2,5,6,11,12,100,110,210,212};
  if (opt_infile == "") {
    const std::string baseLocation = "/mnt/c/Users/sayef/IdeaProjects/traveler-integrated/data_handler/cgal_libs/cgal_server/location_data/";
    opt_infile = baseLocation + "t" + ".loc";
    
  }

  std::ifstream infile(opt_infile);

  if (!infile) {
    std::cerr << "Error: Unable to open file.\n";
    return 1;
  }

  uint64_t value;
  // Read integers line by line
  while (infile >> value) {
      v.push_back(value);
  }

  infile.close();

  
  // npoints = v.size();
  // // computation of condensed distance matrix
  // double* distmat = new double[(npoints*(npoints-1))/2];
  // k = 0;
  // for (i=0; i<npoints; i++) {
  //   for (j=i+1; j<npoints; j++) {
  //     distmat[k] = distance(v[i], v[j]);
  //     k++;
  //   }
  // }

  npoints = v.size() / 2;
  // computation of condensed distance matrix
  double* distmat = new double[(npoints*(npoints-1))/2];
  k = 0;
  for (i=0; i<npoints; i++) {
    for (j=i+1; j<npoints; j++) {
      distmat[k] = distance_event(v[i*2], v[(i*2)+1], v[j*2], v[(j*2)+1]);
      k++;
    }
  }

  // clustering call
  int* merge = new int[2*(npoints-1)];
  double* height = new double[npoints-1];
  int* node_size = new int[2*(npoints-1)];

  double* start_events = new double[2*(npoints-1)];
  double* end_events = new double[2*(npoints-1)];
  hclust_fast(npoints, distmat, opt_method, merge, height, node_size);

  int left_node_index, right_node_index, root_node_index;
  for (i=0; i<npoints-1; i++) {
    left_node_index = changeMerge(merge[i], npoints);
    if(left_node_index < npoints) {
      start_events[left_node_index] = v[left_node_index*2];
      end_events[left_node_index] = v[(left_node_index*2)+1];
    }

    right_node_index = changeMerge(merge[i+npoints-1], npoints);
    if(right_node_index < npoints) {
      start_events[right_node_index] = v[right_node_index*2];
      end_events[right_node_index] = v[(right_node_index*2)+1];
    }

    if(end_events[left_node_index] > start_events[right_node_index]) {
      std::swap(left_node_index, right_node_index);
    }
    
    root_node_index = i + npoints;
    start_events[root_node_index] = start_events[left_node_index];
    end_events[root_node_index] = end_events[right_node_index];
  }

  int* labels = new int[npoints];
  // cutree_k(npoints, merge, 3, labels);
  cutree_cdist(npoints, merge, height, 40000, labels);
  
  // print result
  for (i=0; i<npoints; i++) {
    printf("%d: %.0lf, %.0lf, %d\n",
           i, v[i*2], v[(i*2)+1], labels[i]);
  }
  printf("\n");
  // print result
  for (i=0; i<npoints-1; i++) {

    printf("%.1lf, %.1lf, %0.1lf, %.1lf  --> left(%d %d) right(%d %d) root(%d %d)\n",
           (double)changeMerge(merge[i], npoints), (double)changeMerge(merge[i+npoints-1], npoints), (double)height[i], (double)node_size[i],
           (int)start_events[changeMerge(merge[i], npoints)], (int)end_events[changeMerge(merge[i], npoints)],
           (int)start_events[changeMerge(merge[i+npoints-1], npoints)], (int)end_events[changeMerge(merge[i+npoints-1], npoints)],
           (int)start_events[i+npoints], (int)end_events[i+npoints]);
  }


  std::cout << std::endl;
  // trying the search query
  // search logic
  // if bin size is less than cluster length, then go down
  // else return the start end point of the current cluster
  // dfs on the start and end time query
  int start_time = 10;
  int end_time = 1150;
  int bin_size = 10;
  i = npoints-2;
  dfs(merge, npoints, start_time, end_time, bin_size, npoints-2, start_events, end_events);
  // dfs(merge, npoints, start_time, end_time, bin_size, i+npoints-1, start_events, end_events);

  
  // clean up
  delete[] distmat;
  delete[] merge;
  delete[] height;
  delete[] labels;
  delete[] node_size;
  delete[] start_events;
  delete[] end_events;
  
  return 0;
}

// // main program
// int main_back(int argc, char** argv)
// {

//   int i,j,k,npoints;

//   // parse command line
  // std::string opt_infile;
//   int opt_method = HCLUST_METHOD_SINGLE;
  // const char* usagemsg = "Usage: hclust-demo <infile> [-m (single|complete|average|median)]\n";
  // for (i=1; i<argc; i++) {
  //   if (0 == strcmp(argv[i], "-m")) {
  //     i++;
  //     if (i<argc) {
  //       if (0 == strcmp(argv[i], "single"))
  //           opt_method = HCLUST_METHOD_SINGLE;
  //       else if (0 == strcmp(argv[i], "complete"))
  //           opt_method = HCLUST_METHOD_COMPLETE;
  //       else if (0 == strcmp(argv[i], "average"))
  //           opt_method = HCLUST_METHOD_AVERAGE;
  //       else if (0 == strcmp(argv[i], "median"))
  //           opt_method = HCLUST_METHOD_MEDIAN;
  //       else {
  //         fputs(usagemsg, stderr);
  //         return 1;
  //       }
  //     } else {
  //       fputs(usagemsg, stderr);
  //       return 1;
  //     }
  //   }
  //   else if (argv[i][0] == '-') {
  //     fputs(usagemsg, stderr);
  //     return 1;
  //   }
  //   else {
  //     opt_infile = argv[i];
  //   }
  // }
  // if (opt_infile == "") {
  //   fputs(usagemsg, stderr);
  //   return 1;
  // }
  
//   // read line segments from input file
//   std::vector<Segment> segs;
//   double x1,x2,y1,y2;
//   FILE* f = fopen(opt_infile.c_str(), "r");
//   if (!f) {
//     fprintf(stderr, "Cannot open '%s'\n", opt_infile.c_str());
//     return 2;
//   }
//   npoints = 0;
//   while (!feof(f)) {
//     npoints++;
//     k = fscanf(f, "%lf,%lf,%lf,%lf\n", &x1, &x2, &y1, &y2);
//     if (k != 4) {
//       fprintf(stderr, "Error in line %i of '%s': wrong format\n", npoints, argv[1]);
//       return 3;
//     }
//     segs.push_back(Segment(Point(x1,x2), Point(y1,y2)));
//   }
//   fclose(f);

//   // computation of condensed distance matrix
//   double* distmat = new double[(npoints*(npoints-1))/2];
//   k = 0;
//   for (i=0; i<npoints; i++) {
//     for (j=i+1; j<npoints; j++) {
//       distmat[k] = distance(segs[i], segs[j]);
//       k++;
//     }
//   }

//   // clustering call
//   int* merge = new int[2*(npoints-1)];
//   double* height = new double[npoints-1];
//   hclust_fast(npoints, distmat, opt_method, merge, height);

//   int* labels = new int[npoints];
//   cutree_k(npoints, merge, 5, labels);
//   //cutree_cdist(npoints, merge, height, 0.5, labels);
  
//   // print result
//   for (i=0; i<npoints; i++) {
//     printf("%3.2f,%3.2f,%3.2f,%3.2f,%i\n",
//            segs[i].p1.x, segs[i].p1.y, segs[i].p2.x, segs[i].p2.y, labels[i]);
//   }
  
//   // clean up
//   delete[] distmat;
//   delete[] merge;
//   delete[] height;
//   delete[] labels;

  
//   return 0;
// }
