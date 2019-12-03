#ifndef CLUSTERING_FUNCS
#define CLUSTERING_FUNCS

#include <iostream>
#include <list>
#include <vector>
#include <unordered_map>
#include <random>
#include <chrono>
#include <time.h>
#include <float.h>
#include <algorithm>
#include <random>

#include "my_vector.hpp"
#include "my_curve.hpp"
#include "util.hpp"
#include "lsh.hpp"

template<class T> std::vector<T>* initialization2(std::vector<T>*, unsigned int, double(*distance_metric)(T&, T&));
template<class T> void updatedist(T*, std::list<std::pair<double,T*>> &, double(*distance_metric)(T&, T&));
template<class T> T* newcentroid(std::list<std::pair<double,T*>> &);
template<class T> T* rangebinarysearch(double, std::pair<double,T*>*, int);
template<class T> std::vector<T>* initialization1(std::vector <T>* data, unsigned int k);
template<class T> std::vector<T*>* assigment1(std::vector<T>* data, std::vector<T>* centers,
                                    unsigned int k, double(*distance_metric)(T&, T&));
template<class T> std::vector<T*>* assigment2(std::vector<T>* data, std::vector<T>* centers,
                                    unsigned int k, lsh* lsh_model, double(*distance_metric)(T&, T&));
template<class T> std::vector<T>* update1(std::vector<T>* data, std::vector<T>* centers, std::vector<T*> *clusters,
                                    unsigned int k, lsh* lsh_model, T(*get_mean)(unsigned int, std::vector<T*>&),
                                    double(*distance_metric)(T&, T&));
template<class T> std::vector<T>* update2(std::vector<T>* data, std::vector<T>* centers, std::vector<T*> *clusters,
                                    unsigned int k, T(*get_mean)(unsigned int, std::vector<T*>&));
template<class T> std::vector<T>* update1_brute(std::vector<T>* data, std::vector<T>* centers, std::vector<T*> *clusters,
                                    unsigned int k, double(*distance_metric)(T&, T&));
template<class T> T get_mean(unsigned int dimentions, std::vector<T*> &cluster);
template<class T> double vector_diff(std::vector<T> *vector1, std::vector<T> *vector2);
template<class T> bool old_centers_equal_new_centers(std::vector<T> *old_centers,  std::vector<T> *new_centers, double tolerance);
template<class T> bool old_clusters_equal_new_clusters(std::vector<T*> *old_clusters, std::vector<T*> *new_clusters, unsigned int k);
template<class T> std::list<double> * silhouette(std::vector<T*>* clusters,std::vector<T>* centers, unsigned int k,
                                        unsigned int n, double(*distance_metric)(T&, T&));

#endif
