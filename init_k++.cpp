#include <iostream>
#include <list>
#include <random>
#include <chrono>

#include "my_vector.hpp"
#include "util.hpp"
#include "init_k++.hpp"

using namespace std;

#define DEBUG 0

list <my_vector>* initialization2(list<my_vector> *vectorlist , unsigned int k,double(*distance_metric)(my_vector&, my_vector&)){
  list<my_vector>* centroids = new list<my_vector>;
  pair<double,my_vector*> dummypair;
  list<pair<double,my_vector*>> mindistances;

  //initializing minimum distances array as infinity
  for(list<my_vector>::iterator it=vectorlist->begin(); it!=vectorlist->end(); ++it){
    dummypair.first = numeric_limits<float>::infinity();
    dummypair.second = &(*it);
    mindistances.push_back(dummypair);//creates a copy of the dummypair
  }

  //Choose a vector uniformly at random to be the 1st centroid
  //random generator
  auto seed_t=chrono::system_clock::now().time_since_epoch();
  auto seed_m=chrono::duration_cast<chrono::nanoseconds>(seed_t);
  default_random_engine generator (seed_m.count());
  uniform_int_distribution<int> distribution(0,vectorlist->size()-1);
  //random generator end
  list<my_vector>::iterator rand = vectorlist->begin(); //first point to the 1st element then increment by chance
  advance(rand, distribution(generator)); // now it points to a random element(point)
  centroids->push_back(*rand); //this is the first of the k centroids
  unsigned int t = 1;

  #if DEBUG
    rand->print_vec();
    //std::cout << *(centroids->begin())->get_dimentions() << '\n';
  #endif

  for (; t < k; t++) {//find the rest k-t centroids
    //update the mindistances with the new centroid found previous
    updatedist(&centroids->back(),mindistances, distance_metric);
    centroids->push_back(*newcentroid(mindistances)); //the newcenter (choosen by mindistances)is saved in centroids
  }
  return centroids;
}

// updates the new min_distance with the centroid, ifonlyif it is smaller than the previous min dist
//this function finds for each point the distance from centroid and if it is smallest from the previous min distance
void updatedist(my_vector* centroid,list<pair<double,my_vector*>> &minD,double(*distance_metric)(my_vector&, my_vector&)){
  for(list<pair<double,my_vector*>>::iterator it=minD.begin(); it!=minD.end(); ++it){
    // std::cout << "kkk" << '\n';
    // std::cout << (it->second)->dim << '\n';
    // std::cout << centroid->dim << '\n';
    // std::cout << "kkk" << '\n';
    double new_dist = distance_metric(*centroid,*(it->second));
    if(new_dist< it->first)
      // to remember closest centroid requires tuple insert here
      it->first = new_dist;
    if(it->second==centroid)//double check centroids should not be taken into cosideration//todo
      it->first=0;
  }
}

my_vector* newcentroid(list<pair<double,my_vector*>> &distr){
  pair<double,my_vector*>* prob_array = new pair<double,my_vector*> [distr.size()];//TODO delete this new
  if(prob_array == NULL)
    std::cout << "new didn't allocate the space" << '\n';
  //sum is used to create the probability array
  unsigned int i=0;
  double sum = 0;
  //create the probability array from distribution array
  for(list<pair<double,my_vector*>>::iterator it=distr.begin(); it!=distr.end(); ++it){
    if(it->first != 0){//not a centroid maybe there should be a check in an other way
      sum+=it->first; //todo ^2
      prob_array[i].first = sum;
      prob_array[i].second = it->second;
      i++;
    }
  }
  //fill the rest of the array usually only the centroids with sum+1 could have been infinity
  for (unsigned int j=i; j < distr.size(); j++) {
    prob_array[j].first = sum+1;
    prob_array[j].second = NULL;
  }
  //random generator
  auto seed_t=chrono::system_clock::now().time_since_epoch();
  auto seed_m=chrono::duration_cast<chrono::nanoseconds>(seed_t);
  default_random_engine generator (seed_m.count());
  uniform_real_distribution<double> distribution(0,sum);//how to make it (]
  //select uniformy an element from (0,sum)
  double rng=distribution(generator);
  return rangebinarysearch(rng,prob_array,i);
}

my_vector* rangebinarysearch(double target,pair<double,my_vector*>* p,int r/*it is the length*/){
  int l = 0;
  while (l <= r) {
        int m = l + (r - l) / 2;
        // Check if target is between p[m-1]<target<p[m]
        if (p[m-1].first < target  && target <= p[m].first)
            return p[m].second;
        // If target is greater, ignore left half
        if (p[m].first < target)
            l = m + 1;
        // If target is smaller, ignore right half
        else
            r = m - 1;
    }
  std::cout << "binary search didn't find the target error" << '\n';
  return NULL;
}
