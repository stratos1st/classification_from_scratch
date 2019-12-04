#include "clustering_funcs.hpp"

using namespace std;

#define DEBUG 0

//========================================================== initialization2================================================================
template<class T>
vector <T>* initialization2(vector<T> *vectorlist , unsigned int k, double(*distance_metric)(T&, T&)){
  vector<T>* centroids = new vector<T>;
  pair<double,T*> dummypair;
  list<pair<double,T*>> mindistances;

  //initializing minimum distances array as infinity
  for(auto it=vectorlist->begin(); it!=vectorlist->end(); ++it){
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

  centroids->push_back(vectorlist->at(distribution(generator))); //this is the first of the k centroids (a random one)
  unsigned int t = 1;

  for (; t < k; t++) {//find the rest k-t centroids
    //update the mindistances with the new centroid found previous
    updatedist(&centroids->back(),mindistances, distance_metric);
    centroids->push_back(*newcentroid(mindistances)); //the newcenter (choosen by mindistances)is saved in centroids
  }
  return centroids;
}
template vector<my_vector>* initialization2(vector<my_vector>*, unsigned int, double(*distance_metric)(my_vector&, my_vector&));
template vector<my_curve>* initialization2(vector<my_curve>*, unsigned int, double(*distance_metric)(my_curve&, my_curve&));

// updates the new min_distance with the centroid, ifonlyif it is smaller than the previous min dist
//this function finds for each point the distance from centroid and if it is smallest from the previous min distance
template<class T>
void updatedist(T* centroid, list<pair<double,T*>> &minD, double(*distance_metric)(T&, T&)){
  for(auto it=minD.begin(); it!=minD.end(); ++it){
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
template void updatedist(my_vector*, list<pair<double,my_vector*>> &, double(*distance_metric)(my_vector&, my_vector&));
template void updatedist(my_curve*, list<pair<double,my_curve*>> &, double(*distance_metric)(my_curve&, my_curve&));

template<class T>
T* newcentroid(list<pair<double,T*>> &distr){
  pair<double,T*>* prob_array = new pair<double,T*> [distr.size()];
  //sum is used to create the probability array
  unsigned int i=0;
  double sum = 0;
  //create the probability array from distribution array
  for(auto it=distr.begin(); it!=distr.end(); ++it){
    if(it->first != 0){//not a centroid maybe there should be a check in an other way
      sum+=pow(it->first,2);
      prob_array[i].first = sum;
      prob_array[i].second = it->second;
      i++;
    }
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
template my_vector* newcentroid(list<pair<double,my_vector*>> &);
template my_curve* newcentroid(list<pair<double,my_curve*>> &);

template<class T>
T* rangebinarysearch(double target, pair<double,T*>* p, int r/*it is the length*/){
  int l = 0;
  while (l <= r) {
        int m = l + (r - l) / 2;
        // Check if target is between p[m-1]<target<p[m]
        if (p[m-1].first < target  && target <= p[m].first){
            T* tmp=p[m].second;
            delete[] p;
            return tmp;
        }
        // If target is greater, ignore left half
        if (p[m].first < target)
            l = m + 1;
        // If target is smaller, ignore right half
        else
            r = m - 1;
    }
  std::cout << "binary search didn't find the target error" << '\n';
  exit(1);
}
template my_vector* rangebinarysearch(double, pair<double,my_vector*>*, int);
template my_curve* rangebinarysearch(double, pair<double,my_curve*>*, int);

//========================================================== initialization1================================================================
template<class T>
vector<T>* initialization1(vector <T>* data, unsigned int k){//8eli srand (time(NULL)); apo tin main na kalesti
  // --------------------------------------------------------------------initialization 1 ( nomizo afto lei)
  #if DEBUG
  cout<<"initialization1\n";
  #endif

  vector<int> random_numbers;
  for(unsigned int i=0; i<data->size(); i++)
      random_numbers.push_back(i);
  random_shuffle(random_numbers.begin(), random_numbers.end());

  vector <T>* centers=new vector <T>;
  for(unsigned int i=0;i<k;i++)
    centers->push_back(data->at(random_numbers[i]));

  #if DEBUG
  cout<<"centers\n";
  for(auto i : *centers)
    i.print_vec();
  #endif
  return centers;
}
template vector<my_vector>* initialization1(vector <my_vector>* data, unsigned int k);
template vector<my_curve>* initialization1(vector <my_curve>* data, unsigned int k);

//========================================================== assigment1================================================================
template<class T>
vector<T*>* assigment1(vector<T>* data, vector<T>* centers, unsigned int k, double(*distance_metric)(T&, T&)){
  // --------------------------------------------------------------------asigment 1 (nomizo enoi brute force)
  #if DEBUG
  cout<<"assigment1\n";
  #endif
  vector<T*> *clusters=new vector<T*>[k];
  unsigned int cluster_no=0, best_cluster=0;
  double minn=DBL_MAX,tmp=0.0;
  for(auto &i : *data){
    minn=DBL_MAX;
    cluster_no=0;
    best_cluster=0;
    for(auto &j : *centers){
      cluster_no++;
      tmp=distance_metric(i, j);
      if(tmp<minn){
        best_cluster=cluster_no;
        minn=tmp;
      }
    }
    if(best_cluster==0){
      cerr<<"\n\n!! asigment 1 ERROR cluster not found !!\n\n";
      exit(1);
    }
    clusters[best_cluster-1].push_back(&i);
  }

  #if DEBUG
  //mpori ena cluster na ine keno!!
  for(unsigned int i=0;i<k;i++)
    if(clusters[i].size()==0){
      cout<<"cluster of size 0\n";
      // exit(1);
    }

  cout<<"clusters\n";
  unsigned int ii=0;
  for(unsigned int i=0;i<k;i++){
    cout<<"\tcluster"<<++ii<<endl;
    for(T* j : clusters[i])
      j->print_vec();
  }
  #endif

  return clusters;
}
template vector<my_vector*>* assigment1(vector<my_vector>* data, vector<my_vector>* centers, unsigned int k, double(*distance_metric)(my_vector&, my_vector&));
template vector<my_curve*>* assigment1(vector<my_curve>* data, vector<my_curve>* centers, unsigned int k, double(*distance_metric)(my_curve&, my_curve&));

//========================================================== assigment2================================================================
template<class T>
vector<T*>* assigment2(vector<T>* data, vector<T>* centers, unsigned int k, lsh* lsh_model, double(*distance_metric)(T&, T&)){
// --------------------------------------------------------------------asigment 2
  #if DEBUG
  cout<<"assigment2\n";
  #endif
  vector<T*> *clusters=new vector<T*>[k];
  hash<T*> hasher;


  unordered_map<unsigned int, T*> left_to_classify;
  //perni kapia ora. isos prepi na to ftiaxnoume mia fora stin main kai na pername antigrafo
  //etsi xanoume poli mnimi ala kerdizoume ligo xrono
  for(T& i : *data)
    left_to_classify.insert(make_pair(hasher(&i),&i));

  unordered_map<unsigned int, T*> *close_points[k];
  for(unsigned int i=0;i<k;i++){
    close_points[i]=lsh_model->find_bucket(centers->at(i), distance_metric);
    #if DEBUG
    cout<<"lsh buckets\n-----------"<<i<<endl;
    for (auto v : *close_points[i])
      v.second->print_vec();
    #endif
  }

  /*create uper triangular table with no diagonal
    the table is filled like the example
    ex k=5
      0 1 2 3 4
     -----------
    0|
    1|.
    2|. .
    3|. . .
    4|. . . .

  */
  unordered_map<unsigned int, T*> *intersection[k-1];
  for(unsigned int i=0;i<k-1;i++)
    intersection[i]= new unordered_map<unsigned int, T*>[i+1];
  //find intersection of close_points
  unordered_map<unsigned int, T*> *big, *small;
  for(unsigned int ii=0;ii<k-1;ii++)
    for(unsigned int j=0;j<=ii;j++){//for every . in intersection
      big=close_points[ii+1];
      small=close_points[j];
      if(big->size()<small->size())
        swap(big,small);//iterate throw the smallest set
      for(auto i = small->begin(); i != small->end(); i++)//for every element in the set
        if(big->find(hasher(i->second)) != big->end())//if the element exists in the other set
          intersection[ii][j].insert(*i);//insert it in intersection
    }

  #if DEBUG
  cout<<"intersection\n";
  for(unsigned int i=0;i<k-1;i++)
    for(unsigned int j=0;j<=i;j++){
      cout<<"\t"<<i+1<<" "<<j<<endl;
      for (auto v : intersection[i][j])
            v.second->print_vec();
    }
  #endif

  double min_dist,tmp1,tmp2;
  unsigned int closest_center;

  //classify all intersection elements
  #if DEBUG
  cout<<"-all intersection elements\n";
  #endif
  for(unsigned int i=0;i<k-1;i++){
    for(unsigned int j=0;j<=i;j++){
      for(auto v=intersection[i][j].begin();v!=intersection[i][j].end();++v){//for every element in intersection
        if(left_to_classify.find(hasher(v->second)) != left_to_classify.end()){//if the element has yet to be classified
          min_dist=DBL_MAX;
          for(unsigned int ik=0;ik<k-1;ik++){
            for(unsigned int jk=0;jk<=ik;jk++){//search all the intersections
              if(intersection[ik][jk].find(hasher(v->second)) !=  intersection[ik][jk].end()){//if it exists in any of them find the min dist and the closest_center
                tmp1=distance_metric(centers->at(ik+1), *v->second);
                tmp2=distance_metric(centers->at(jk), *v->second);
                if(tmp1<min_dist){
                  min_dist=tmp1;
                  closest_center=ik+1;
                }
                if(tmp2<min_dist){
                  min_dist=tmp2;
                  closest_center=jk;
                }
              }
            }
          }
          #if DEBUG
          v->second->print_vec();
          #endif
          clusters[closest_center].push_back(v->second);//push the element in the closest cluster
          left_to_classify.erase(hasher(v->second));//erase it from left_to_classify
        }
      }
    }
  }


  //classify all close_points elements (exept the ones classified in intersection)
  #if DEBUG
  cout<<"-all close_points elements (exept the ones classified in intersection)\n";
  #endif
  for(unsigned int i=0;i<k;i++)
    for (auto v=close_points[i]->begin();v!=close_points[i]->end();++v)//for all the close points
      if(left_to_classify.find(hasher(v->second)) != left_to_classify.end()){//if the element has yet to be classified
        #if DEBUG
        v->second->print_vec();
        #endif
        clusters[i].push_back(v->second);//push the element in the closest cluster
        left_to_classify.erase(hasher(v->second));//erase it from left_to_classify
      }

  //classify all the rest
  #if DEBUG
  cout<<"-all the rest\n";
  #endif
  for (auto v=left_to_classify.begin();v!=left_to_classify.end();){//for every element in left_to_classify
    min_dist=DBL_MAX;
    for(unsigned int i=0;i<k;i++){//check all the centers
      tmp1=distance_metric(*v->second, centers->at(i));
      if(min_dist>tmp1){//find the cosest
        min_dist=tmp1;
        closest_center=i;
      }
    }
    clusters[closest_center].push_back(v->second);//push the element in the closest cluster
    #if DEBUG
    v->second->print_vec();
    v=left_to_classify.erase(v);
    #else
    ++v;
    #endif

  }

  #if DEBUG
  if(left_to_classify.size()!=0){
    cerr<<"\n\nleft_to_classify.size()!=0\n\n";
    for (auto v : left_to_classify)
      v.second->print_vec();
    exit(1);
  }

  //mpori ena cluster na ine keno!!
  for(unsigned int i=0;i<k;i++)
    if(clusters[i].size()==0){
      cout<<"cluster of size 0\n";
      // exit(1);
    }

  cout<<"clusters\n";
  unsigned int ii=0;
  for(unsigned int i=0;i<k;i++){
    cout<<"\tcluster"<<++ii<<endl;
    for(T* j : clusters[i])
      j->print_vec();
  }
  #endif

  for(unsigned int i=0;i<k-1;i++){
    intersection[i]->clear();
    delete[] intersection[i];
  }
  for(unsigned int i=0;i<k;i++){
    close_points[i]->clear();
    delete close_points[i];
  }

  return clusters;
}
template vector<my_vector*>* assigment2(vector<my_vector>* data, vector<my_vector>* centers, unsigned int k, lsh* lsh_model, double(*distance_metric)(my_vector&, my_vector&));
template vector<my_curve*>* assigment2(vector<my_curve>* data, vector<my_curve>* centers, unsigned int k, lsh* lsh_model, double(*distance_metric)(my_curve&, my_curve&));

//========================================================== update1================================================================
template<class T>
vector<T>* update1(vector<T>* data, vector<T>* centers, vector<T*> *clusters, unsigned int k, lsh* lsh_model,
                           T(*get_mean)(unsigned int, vector<T*>&), double(*distance_metric)(T&, T&)){//nomizo oti ine afto. afto ine genika kalo ke to kratame gt maresi san bonus
  // --------------------------------------------------------------------update 1
  #if DEBUG
  cout<<"update1\n";
  #endif
  vector<T>* new_centers=new vector<T>;

  for(unsigned int i=0;i<k;i++)
    if(clusters[i].size()!=0){
      T average(get_mean(data->front().get_dimentions(),clusters[i]));
      pair<T*,double> tmp=lsh_model->find_NN(average, distance_metric);
      if(tmp.second==DBL_MAX){
        #if DEBUG
        cerr<<"wrong lsh arguments!! no NN found!! going brute\n";
        #endif
        vector<T> *center_tmp=new vector<T>;
        center_tmp->push_back(centers->at(i));
        vector<T>* tmp_vec=update1_brute(data,center_tmp,&clusters[i],1,distance_metric);
        new_centers->push_back(tmp_vec->at(0));
        tmp_vec->clear();
        delete tmp_vec;
        center_tmp->clear();
        delete center_tmp;
      }
      else
        new_centers->push_back(*tmp.first);
    }
    else
      new_centers->push_back(centers->at(i));

  #if DEBUG
  cout<<"updated centers\n";
  for(auto i : *new_centers)
    i.print_vec();
  #endif

  return new_centers;
}
template vector<my_vector>* update1(vector<my_vector>* data, vector<my_vector>* centers, vector<my_vector*> *clusters, unsigned int k, lsh* lsh_model,
                          my_vector(*get_mean)(unsigned int, vector<my_vector*>&), double(*distance_metric)(my_vector&, my_vector&));
template vector<my_curve>* update1(vector<my_curve>* data, vector<my_curve>* centers, vector<my_curve*> *clusters, unsigned int k, lsh* lsh_model,
                          my_curve(*get_mean)(unsigned int, vector<my_curve*>&), double(*distance_metric)(my_curve&, my_curve&));

//========================================================== update1_brute================================================================
template<class T>
vector<T>* update1_brute(vector<T>* data, vector<T>* centers, vector<T*> *clusters, unsigned int k, double(*distance_metric)(T&, T&)){
  // --------------------------------------------------------------------update 1
  #if DEBUG
  cout<<"update1_brute\n";
  #endif
  vector<T>* new_centers=new vector<T>;

  double minn,sum;
  T* best_center;
  for(unsigned int i=0;i<k;i++)
    if(clusters[i].size()!=0){
      minn=DBL_MAX;
      for(auto j: clusters[i]){
        sum=0;
        for(auto jj: clusters[i])
          sum+=distance_metric(*j,*jj);
        if(minn>sum){
          minn=sum;
          best_center=j;
        }
      }
      new_centers->push_back(*best_center);
    }
    else
      new_centers->push_back(centers->at(i));

  #if DEBUG
  cout<<"updated centers\n";
  for(auto i : *new_centers)
    i.print_vec();
  #endif

  return new_centers;
}
template vector<my_vector>* update1_brute(vector<my_vector>* data, vector<my_vector>* centers, vector<my_vector*> *clusters, unsigned int k, double(*distance_metric)(my_vector&, my_vector&));
template vector<my_curve>* update1_brute(vector<my_curve>* data, vector<my_curve>* centers, vector<my_curve*> *clusters, unsigned int k, double(*distance_metric)(my_curve&, my_curve&));

//========================================================== update2================================================================
template<class T>
vector<T>* update2(vector<T>* data, vector<T>* centers, vector<T*> *clusters, unsigned int k, T(*get_mean)(unsigned int, vector<T*>&)){
  // --------------------------------------------------------------------update 2 (nomizo enoi brute force)
  vector<T>* new_centers=new vector<T>;

  for(unsigned int i=0;i<k;i++)
    if(clusters[i].size()!=0)
      new_centers->push_back(get_mean(data->front().get_dimentions(),clusters[i]));
    else
      new_centers->push_back(centers->at(i));

  #if DEBUG
  cout<<"updated centers\n";
  for(auto i : *new_centers)
    i.print_vec();
  #endif

  return new_centers;
}
template vector<my_vector>* update2(vector<my_vector>* data, vector<my_vector>* centers, vector<my_vector*> *clusters, unsigned int k, my_vector(*get_mean)(unsigned int, vector<my_vector*>&));
template vector<my_curve>* update2(vector<my_curve>* data, vector<my_curve>* centers, vector<my_curve*> *clusters, unsigned int k, my_curve(*get_mean)(unsigned int, vector<my_curve*>&));

//========================================================== exit conditions================================================================
template<class T>
bool old_clusters_equal_new_clusters(vector<T*> *old_clusters, vector<T*> *new_clusters, unsigned int k){
  for(unsigned int i=0;i<k;i++)
    if(old_clusters[i]!=new_clusters[i])
      return false;
  return true;
}
template bool old_clusters_equal_new_clusters<>(vector<my_vector*> *old_clusters, vector<my_vector*> *new_clusters, unsigned int k);
template bool old_clusters_equal_new_clusters<>(vector<my_curve*> *old_clusters, vector<my_curve*> *new_clusters, unsigned int k);

template<class T>
bool old_centers_equal_new_centers(vector<T> *old_centers,  vector<T> *new_centers, double tolerance){
  return vector_diff(old_centers,new_centers)<tolerance;
}
template bool old_centers_equal_new_centers<>(vector<my_vector> *old_centers,  vector<my_vector> *new_centers, double tolerance);
template bool old_centers_equal_new_centers<>(vector<my_curve> *old_centers,  vector<my_curve> *new_centers, double tolerance);

//========================================================== mean functions================================================================
template<> my_vector get_mean<my_vector>(unsigned int dimentions, vector<my_vector*> &cluster){//DEN prepi to cluster na ine 0!!
  #if DEBUG
  cout<<"get_vector_mean\n";
  #endif
  my_vector average(dimentions);

  for(unsigned int i=0;i<dimentions;i++)
  average.coordinates[i]=0;

  for(my_vector* j : cluster)
  for(unsigned int i=0;i<dimentions;i++)
  average.coordinates[i]+=j->coordinates[i];
  //DEN mpori na ine 0 to clusters[i].size()
  for(unsigned int i=0;i<dimentions;i++)
  average.coordinates[i]/=cluster.size();

  return average;
}

template<> my_curve get_mean<my_curve>(unsigned int dimentions, vector<my_curve*> &cluster){//DEN prepi to cluster na ine 0!!
  #if DEBUG
  cout<<"get_curve_mean\n";
  #endif
  //bres to l
  double average=0.0;
  for(my_curve* i: cluster)
    average+=i->numofvectors;
  average/=cluster.size();//DEN mpori na ine 0!!
  unsigned int l=(unsigned int)floor(average);
  //arxikopiise to C
  my_curve *c;
  for(my_curve* it: cluster)//FIXME perno to proto me length==l eno prepi na perno tixeo
    if(it->numofvectors>=average){
      c=new my_curve(l,dimentions);
      unsigned int jj=rand()%(it->numofvectors-l+1);//ksekina apo mia 8esi
      for(unsigned int j=0;j<l;j++)//ke pigene mexri jj+l
        for(unsigned int i=0;i<dimentions;i++)//gia ka8e diastasi
          c->vectors[j]->coordinates[i]=it->vectors[jj+j]->coordinates[i];//antegrapse sto C
      break;
    }

  my_curve c2(l,dimentions);
  vector<my_vector*>* a=new vector<my_vector*>[l];
  list<pair<unsigned int,unsigned int>>* ipairs;
  while(1){
    c2=*c;
    for(unsigned int i=0;i<l;i++)
      a[i].clear();
    for(my_curve* it: cluster){
      ipairs=MinMatching(*c,*it);
      for(pair<unsigned int,unsigned int> ij: *ipairs)
        a[ij.first-1].push_back(it->vectors[ij.second-1]);
      ipairs->clear();
      delete ipairs;
    }
    for(unsigned int i=0;i<l;i++)
      *c->vectors[i]=get_mean(dimentions,a[i]);
    if(c2==*c)
      break;
    //delete c2;
  }

  //delete c2;
  for(unsigned int i=0;i<l;i++)
    a[i].clear();
  delete[] a;

  my_curve ans(*c);
  delete c;
  return ans;
}

//========================================================== helping functions================================================================
template<> double vector_diff<my_vector>(vector<my_vector> *list1, vector<my_vector> *list2){
  double maxx=0;
  vector<my_vector>::iterator i=list1->begin();

  if(list1->size()!=list2->size()){
    cerr<<"\n\n!!error list_diff vectors !=size !!\n\n";
    exit(1);
  }

  for(my_vector j : *list2){
    maxx=max(manhattan_distance(*i,j),maxx);
    i++;
  }
  return maxx;
}

template<> double vector_diff<my_curve>(vector<my_curve> *list1, vector<my_curve> *list2){
  double maxx=0;
  vector<my_curve>::iterator i=list1->begin();

  if(list1->size()!=list2->size()){
    cerr<<"\n\n!!error list_diff curves !=size !!\n\n";
    exit(1);
  }

  for(my_curve& j : *list2){
    maxx=max(Dtw(*i,j),maxx);
    ++i;
  }
  return maxx;
}

//========================================================== silhouette================================================================
template<class T>
list<double> * silhouette(vector<T*>* clusters,vector<T>* centers,unsigned int k,unsigned int n,double(*distance_metric)(T&, T&)){
  double a[n];
  double b[n];
  unsigned int clustersize = 0;
  unsigned int i = 0;//this is for a[i]
  unsigned int r = 0;//this is for b[r]
  double clustercoef[k];
  double **distarray;
  for (unsigned int x  = 0; x < n; x++) {
    b[x] = 0.0;
    a[x] = 0.0;
  }
  for (unsigned int x  = 0; x < k; x++)
    clustercoef[x] = 0.0;

  //NearestC & 2ndNearestC
  for(unsigned int j=0; j<k; j++){//take a center and its cluster(they sould be sync)

    clustersize = clusters[j].size();
    #if DEBUG
      std::cout << "j="<<j<<"k="<<k <<"i,r = "<<i<<","<<r <<"size="<<clustersize<<"with center"<<(*centers)[j].id<< '\n';
    #endif

    //create lower triangular matrix
    if(clustersize<=1)
      continue;
    distarray = new double* [clustersize-1];//create distance array:for cluster
    for (unsigned int v = 0; v<clustersize-1; v++) {//z:0->clustersize-2 , lenght clustersize-1
      distarray[v]=new double[v+1];
      for (unsigned int l = 0; l < v+1; l++) {//fill all z blocks of mem with doubles
        distarray[v][l] = distance_metric(*clusters[j][v],*clusters[j][l+1]);
      }
    }//done with assigment and creation of dist array
    //---------------------------------------------------compute a[i] for each i in j
    for (unsigned int l = 0; l < clustersize; l++) {
      for (unsigned int m= 0; m < clustersize; m++) {
        if(l>m)
          a[i] += distarray[l-1][m];
        else if(l<m)
          a[i] += distarray[m-1][l];
      }
      if (clustersize>1)
        a[i]=a[i]/(double)(clustersize-1);
      else
        a[i]=0;
      i++;
    }
    //delete dist array
    for (size_t l = 0; l < clustersize-1; l++)
      delete[] distarray[l];
    delete[] distarray;
    //--------------compute b[i] for each i in j-------------------------------
    for (unsigned int v = 0; v < clustersize; v++) {//take a vector of j cluster
      unsigned int secondcenter = k ;// initailized out of bounds on perpuse
      double min=DBL_MAX;
      for (unsigned int c = 0; c < k; c++) {// take a center
        if (c!=j ) {//exclude the center we are in  //maybe add z!=j
          if(distance_metric(*clusters[j][v],(*centers)[c])<min){//find the second best for that v
            min = distance_metric(*clusters[j][v],(*centers)[c]);
            secondcenter = c;
            #if DEBUG
              std::cout << "found the secondcenter of "<<clusters[j][v]->id<<" == "<<(*centers)[c].id<< '\n';
            #endif
          }
        }
      }
      //found second center
      if(secondcenter < k ){//find b[r]
        unsigned int secondcentersize = clusters[secondcenter].size();
        for (unsigned int v2= 0; v2 < secondcentersize; v2++) {//take a vector from the 2nd best cluster
          b[r]+=distance_metric(*clusters[j][v],*clusters[secondcenter][v2]);//distance from v to every other v on the second best cluster
        }
        b[r]=b[r]/(double)secondcentersize;
        r++;
      }
      else{
        std::cout << "!!error did not find 2nd center of vector "<<clusters[j][v]->id << '\n';
        exit(1);
      }
    }
  }

  // ---------------------------------------cumpute clustercoef[j]
  clustersize = 0;
  i = 0;   // i goes from 0->n because sum(clustersize[j]) over j:1->k == n
  unsigned int sum = 0;
  for (unsigned int j = 0; j < k; j++) {
    clustersize = clusters[j].size();
    if(clustersize==0){
      clustercoef[j]=-2;
      continue;
    }
    else if(clustersize==1){
      clustercoef[j]=0;
      continue;
    }
    for (; i < sum+clustersize; i++) {
      // this is s(i)
      if(a[i]<b[i]){
        #if DEBUG
        std::cout << b[i]<<"+"<<a[i] << '\n';
        std::cout << 1-(a[i]/b[i]) << '\n';
        #endif
        clustercoef[j]+=1-(a[i]/b[i]);
      }
      else if(a[i]>b[i])
      {
        #if DEBUG
        std::cout << b[i]<<"-"<<a[i] << '\n';
        std::cout << (b[i]/a[i])-1 << '\n';
        #endif
        clustercoef[j]+=(b[i]/a[i])-1;
      }
    }
    sum+=clustersize;
    clustercoef[j]/=(double)(clustersize);
  }
  #if DEBUG
    for (unsigned int x = 0; x < n; x++) {
      cout<<x+1<<". "<<b[x]<<endl;
    }
    for (unsigned int x = 0; x < n; x++) {
      cout<<x+1<<". "<<a[x]<<endl;
    }
    std::cout << "\nCluster Coefficients" << '\n';
    std::cout << "--------------------" << '\n';
    for (unsigned int x = 0; x < k; x++) {
      cout<<x+1<<". "<<clustercoef[x]<<endl;
    }
  #endif

  list<double> *result = new list<double>;
  for (unsigned int i = 0; i < k; i++)
    result->push_back(clustercoef[i]);

  return result;
}
template list<double> * silhouette(vector<my_vector*>* clusters,vector<my_vector>* centers,unsigned int k,unsigned int n,double(*distance_metric)(my_vector&, my_vector&)=manhattan_distance);
template list<double> * silhouette(vector<my_curve*>* clusters,vector<my_curve>* centers,unsigned int k,unsigned int n,double(*distance_metric)(my_curve&, my_curve&));
