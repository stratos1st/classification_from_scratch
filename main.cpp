#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <float.h>
#include <algorithm>
#include <random>
#include <unistd.h>
#include <string.h>

#include "util.hpp"
#include "my_vector.hpp"
#include "my_curve.hpp"
#include "lsh.hpp"
#include "init_k++.hpp"

using namespace std;

#define DEBUG 0

template<class T>
vector<T>* initialization1(vector <T>* data, unsigned int k);
template vector<my_vector>* initialization1<>(vector <my_vector>* data, unsigned int k);
template vector<my_curve>* initialization1<>(vector <my_curve>* data, unsigned int k);

template<class T>
vector<T*>* assigment1(vector<T>* data, vector<T>* centers, unsigned int k, double(*distance_metric)(T&, T&));
template vector<my_vector*>* assigment1<>(vector<my_vector>* data, vector<my_vector>* centers, unsigned int k, double(*distance_metric)(my_vector&, my_vector&));
template vector<my_curve*>* assigment1<>(vector<my_curve>* data, vector<my_curve>* centers, unsigned int k, double(*distance_metric)(my_curve&, my_curve&));

template<class T>
vector<T*>* assigment2(vector<T>* data, vector<T>* centers, unsigned int k, lsh* lsh_model, double(*distance_metric)(T&, T&));
template vector<my_vector*>* assigment2<>(vector<my_vector>* data, vector<my_vector>* centers, unsigned int k, lsh* lsh_model, double(*distance_metric)(my_vector&, my_vector&));
template vector<my_curve*>* assigment2<>(vector<my_curve>* data, vector<my_curve>* centers, unsigned int k, lsh* lsh_model, double(*distance_metric)(my_curve&, my_curve&));

template<class T>
vector<T>* update1(vector<T>* data, vector<T>* centers, vector<T*> *clusters, unsigned int k, lsh* lsh_model,
                           T(*get_mean)(unsigned int, vector<T*>&), double(*distance_metric)(T&, T&));
template vector<my_vector>* update1<>(vector<my_vector>* data, vector<my_vector>* centers, vector<my_vector*> *clusters, unsigned int k, lsh* lsh_model,
                          my_vector(*get_mean)(unsigned int, vector<my_vector*>&), double(*distance_metric)(my_vector&, my_vector&));
template vector<my_curve>* update1<>(vector<my_curve>* data, vector<my_curve>* centers, vector<my_curve*> *clusters, unsigned int k, lsh* lsh_model,
                          my_curve(*get_mean)(unsigned int, vector<my_curve*>&), double(*distance_metric)(my_curve&, my_curve&));

template<class T>
vector<T>* update2(vector<T>* data, vector<T>* centers, vector<T*> *clusters, unsigned int k, T(*get_mean)(unsigned int, vector<T*>&));
template vector<my_vector>* update2<>(vector<my_vector>* data, vector<my_vector>* centers, vector<my_vector*> *clusters, unsigned int k, my_vector(*get_mean)(unsigned int, vector<my_vector*>&));
template vector<my_curve>* update2<>(vector<my_curve>* data, vector<my_curve>* centers, vector<my_curve*> *clusters, unsigned int k, my_curve(*get_mean)(unsigned int, vector<my_curve*>&));

my_curve get_curve_mean(unsigned int dimentions, vector<my_curve*> &cluster);
my_vector get_vector_mean(unsigned int dimentions, vector<my_vector*> &cluster);

double vector_diff(vector<my_vector> *vector1, vector<my_vector> *vector2);
double vector_diff(vector<my_curve> *vector1, vector<my_curve> *vector2);

template<class T>
bool old_centers_equal_new_centers(vector<T> *old_centers,  vector<T> *new_centers, double tolerance);
template bool old_centers_equal_new_centers<>(vector<my_vector> *old_centers,  vector<my_vector> *new_centers, double tolerance);
template bool old_centers_equal_new_centers<>(vector<my_curve> *old_centers,  vector<my_curve> *new_centers, double tolerance);

template<class T>
bool old_clusters_equal_new_clusters(vector<T*> *old_clusters,  vector<T*> *new_clusters, unsigned int k);
template bool old_clusters_equal_new_clusters<>(vector<my_vector*> *old_clusters,  vector<my_vector*> *new_clusters, unsigned int k);
template bool old_clusters_equal_new_clusters<>(vector<my_curve*> *old_clusters,  vector<my_curve*> *new_clusters, unsigned int k);


//filippos was here
template<class T>
list<double> * vectorizeandsilhouette(list<T*>* lclusters,list<T>* lcenters,unsigned int k,unsigned int n,double(*distance_metric)(T&, T&));
template  list<double> * vectorizeandsilhouette(list<my_vector*>* lclusters,list<my_vector>* lcenters,unsigned int k,unsigned int n,double(*distance_metric)(my_vector&, my_vector&)=manhattan_distance);
template  list<double> * vectorizeandsilhouette(list<my_curve*>* lclusters,list<my_curve>* lcenters,unsigned int k,unsigned int n,double(*distance_metric)(my_curve&, my_curve&));
template<class T>
list<double> * silhouette(vector<T*>* clusters,vector<T>* centers,unsigned int k,unsigned int n,double(*distance_metric)(T&, T&));
template list<double> * silhouette(vector<my_vector*>* clusters,vector<my_vector>* centers,unsigned int k,unsigned int n,double(*distance_metric)(my_vector&, my_vector&)=manhattan_distance);
template list<double> * silhouette(vector<my_curve*>* clusters,vector<my_curve>* centers,unsigned int k,unsigned int n,double(*distance_metric)(my_curve&, my_curve&));

int main(int argc, char** argv){
  srand (time(NULL));
  char input_file[100]("./Input/DataVectors_5_500x100.csv"),out_file[100]("my.out"),options_file[100]("cluster.conf");
  unsigned int k=4, max_iterations=50,lsh_window=6000,g_no=4,grids_no=4,container_sz=10,lsh_l=4,max_curve_sz=10;
  bool vector_input=true,stop_when_centers=true;;
  short int function_matrix[3],complete_flag=0;
  double center_tol=0.01,pad_number=99999.99999;
  my_curve::curve_tol=0.01;
  my_vector::vector_tol=0.01;

  cout<<"arxisame!!\n";

  //------------------------------------parse arguments
  int opt;
  while((opt = getopt(argc, argv, "i:o:c:a:"))!=-1){
    switch(opt){
      case 'i':
        strcpy(input_file,optarg);
        break;
      case 'o':
        strcpy(out_file,optarg);
        break;
      case 'c':
        strcpy(options_file,optarg);
        break;
      case 'a':
        complete_flag=atoi(optarg);//FIXME not working
        break;
      default:
        cout<<"!! WRONG ARGUMENTS !!\n";
        exit(1);
    }
  }

//-----------------------------------------------------read cluster.conf
  ifstream infile(options_file);
  char str[100];
  if (infile.good()){
    while(infile.getline(str,100)){
      sscanf(str,"number_of_clusters: %u",&k);//FIXME does not take k
      sscanf(str,"number_of_grids: %u",&grids_no);
      sscanf(str,"number_of_vector_hash_tables: %u",&lsh_l);
      sscanf(str,"number_of_vector_hash_functions: %u",&g_no);
      sscanf(str,"lsh_window: %u",&lsh_window);
      sscanf(str,"input_contains_vectors: %u",&vector_input);
      sscanf(str,"lsh_multimap_container_size: %u",&container_sz);
      sscanf(str,"stop_when_centers_dont_change: %u",&stop_when_centers);
      sscanf(str,"curve_tolerance: %lf",&my_curve::curve_tol);
      sscanf(str,"vector_tolerance: %lf",&my_vector::vector_tol);
      sscanf(str,"center_tolerance: %lf",&center_tol);
      sscanf(str,"max_iterations: %u",&max_iterations);
      sscanf(str,"lsh_curves_pad_number: %lf",&pad_number);
      sscanf(str,"max_curve_size: %u",&max_curve_sz);
    }
  }
  else{
    cerr << "\n\n!! .conf FILE ERROR !!\n\n";
    exit(1);
  }
  infile.close();

  k=4;
  //cout parameters
  cout<<"starting parameters:"<<"\n\tnumber_of_clusters= "<<k<<"\n\tnumber_of_grids= "<<grids_no
  <<"\n\tlsh_l= "<<lsh_l<<"\n\tnumber_of_vector_hash_functions= "<<g_no<<"\n\tlsh_window= "<<lsh_window
  <<"\n\tvector_input= "<<vector_input<<"\n\tcontainer_sz= "<<container_sz<<"\n\tpad_number= "<<pad_number
  <<"\n\tmax_curve_size= "<<max_curve_sz
  <<"\n\tstop_when_centers_dont_change= "<<stop_when_centers<<"\n\tcurve_tolerance= "<<my_curve::curve_tol
  <<"\n\tvector_tolerance= "<<my_vector::vector_tol<<"\n\tcenter_tolerance= "<<center_tol
  <<"\n\tmax_iterations= "<<max_iterations<<"\n\tinput_file= "<<input_file
  <<"\n\tout_file= "<<out_file<<"\n\tcomplete_flag= "<<complete_flag<<endl;

////-----------------------------------------------------read input file
  // list <my_vector>* data_tmp=read_vector_file(input_file);
  list <my_curve>* data_tmp=read_curve_file(input_file);

  cout<<"reading done\n";
  #if DEBUG
  cout<<"data\n";
  for(auto i : *data_tmp)
    i.print_vec();
  #endif

  lsh *lsh_model;
  double (*distance_metric)(my_curve&,my_curve&);
  my_curve (*mean_function)(unsigned int, vector<my_curve*> &);
  if(vector_input){
    lsh_model=new lsh_vector(data_tmp->front().get_dimentions(),lsh_l,lsh_window,grids_no,container_sz);
    //distance_metric=manhattan_distance;
    // mean_function=get_vector_mean;
  }
  else{
    lsh_model=new lsh_curve(data_tmp->front().get_dimentions(),max_curve_sz,lsh_l,lsh_window,grids_no,pad_number,container_sz);
    distance_metric=Dtw;
    mean_function=get_curve_mean;
  }

  lsh_model->train(data_tmp);
  cout<<"lsh training done\n";

  vector<my_curve>* data=new vector<my_curve>;
  for(auto i: *data_tmp)
    data->push_back(i);
  data_tmp->clear();
  delete data_tmp;

  vector<my_curve> *old_centers,*centers;
  vector<my_curve*> *old_clusters,*clusters;
  unsigned int max_iterations_const=max_iterations;

  for(function_matrix[0]=0;function_matrix[0]<=1;function_matrix[0]++){
    for(function_matrix[1]=0;function_matrix[1]<=1;function_matrix[1]++){
      for(function_matrix[2]=0;function_matrix[2]<=1;function_matrix[2]++){

        //cout algorithms used
        cout<<endl;
        if(function_matrix[0])
          cout<<"initialization1 + ";
        else
          cout<<"initialization2 + ";
        if(function_matrix[1])
          cout<<"assigment1 + ";
        else
          cout<<"assigment2 + ";
        if(function_matrix[2])
          cout<<"update1\n";
        else
          cout<<"update2\n";
        cout<<endl;

        //-----------------------------------------------------------------------initialization
        if(function_matrix[0])
          centers=initialization1(data,k);
        else
          centers=initialization1(data,k);

        //-----------------------------------------------------------------------assigment
        if(function_matrix[1])
          clusters=assigment1(data,centers,k,distance_metric);
        else
          clusters=assigment2(data,centers,k,lsh_model,distance_metric);
        //me tis epiloges gia initialization pou dini ipoti8ete oti ta cluster 8a exoun toulaxiston 1 simio
        //kanonika 8a prepi na to ele3oume prin kalesoume assigment

        unsigned int iteration=0;
        max_iterations=max_iterations_const;
        while(max_iterations--!=0){
          #if DEBUG
          cout<<"\niteration="<<iteration<<endl;
          #endif

          old_centers=centers;
          //----------------------------------------------------------------------update
          if(function_matrix[2])
            centers=update1(data, centers, clusters,k,lsh_model,mean_function,distance_metric);
          else
            centers=update2(data, centers, clusters,k,mean_function);

          // sin8iki break:: nea kentra konta sta palia
          if(stop_when_centers)
            if(old_centers_equal_new_centers(old_centers,centers,center_tol))
              break;
          //delete old_centers
          old_centers->clear();
          delete old_centers;

          old_clusters=clusters;
          //---------------------------------------------------------------------assigment
          if(function_matrix[2])
            clusters=assigment1(data,centers,k,distance_metric);
          else
            clusters=assigment2(data,centers,k,lsh_model,distance_metric);

          //sin8iki break:: nea cluster idia me palia
          if(!stop_when_centers)
            if(old_clusters_equal_new_clusters(old_clusters,clusters,k))
              break;
          //delete old_clusters
          for(unsigned int i=0;i<k;i++)
            old_clusters[i].clear();
          delete[] old_clusters;

          //sin8iki break:: max iterations reached
          if(max_iterations==1){
            cout<<"max iterations reached\n";
            break;
          }
          iteration++;
        }

        //-----------------------------------------------------------------------success
        cout<<"success!!\niteration="<<++iteration<<endl;
        cout<<"end centers\n";
        for(auto i : *centers)
          // i.print_vec();
          cout<<i.id<<endl;
        cout<<"clusters\n";
        unsigned int ii=0;
        for(unsigned int i=0;i<k;i++){
          cout<<"\tcluster"<<++ii<<endl;
          if(clusters[i].size()==0)
            cout<<"cluster of size 0\n";
          // for(auto* j : clusters[i])
          //   j->print_vec();
        }
        //TODO 987897678y7Q($%$*^&%$^#%#@%$#^$&%$)
        // list<double> *a = vectorizeandsilhouette(clusters,centers,centers->size(),data->size(),manhattan_distance);
        // a->clear();
        // delete a;


        //deletes
        if(max_iterations!=1){
          if(stop_when_centers){
            old_centers->clear();
            delete old_centers;
          }
          else{
            for(unsigned int i=0;i<k;i++)
              old_clusters[i].clear();
            delete[] old_clusters;
          }
        }
        centers->clear();
        delete centers;
        for(unsigned int i=0;i<k;i++)
          clusters[i].clear();
        delete[] clusters;

      }
    }
  }

  delete lsh_model;
  data->clear();
  delete data;


  cout<<"\n\nEND!!\n\n";

  return 0;
}

template<class T>
bool old_clusters_equal_new_clusters(vector<T*> *old_clusters, vector<T*> *new_clusters, unsigned int k){
  for(unsigned int i=0;i<k;i++)
    if(old_clusters[i]!=new_clusters[i])
      return false;
  return true;
}

template<class T>
bool old_centers_equal_new_centers(vector<T> *old_centers,  vector<T> *new_centers, double tolerance){
  return list_diff(old_centers,new_centers)<tolerance;
}

template<class T>
vector<T>* initialization1(vector <T>* data, unsigned int k){//8eli srand (time(NULL)); apo tin main na kalesti
  // --------------------------------------------------------------------initialization 1 ( nomizo afto lei)
  #if DEBUG
  cout<<"initialization1\n";
  #endif

  //TODO mpori na iparxi kalieros tropos
  vector<int> random_numbers;
  for(unsigned int i=0; i<data->size(); i++)
      random_numbers.push_back(i);
  random_shuffle(random_numbers.begin(), random_numbers.end());

  vector <T>* centers=new vector <T>;
  for(unsigned int i=0;i<k;i++)
    centers->push_back(data->at(random_numbers[i]%data->size()));

  #if DEBUG
  cout<<"centers\n";
  for(auto i : *centers)
    i.print_vec();
  #endif
  return centers;
}

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

template<class T>
vector<T*>* assigment2(vector<T>* data, vector<T>* centers, unsigned int k, lsh* lsh_model, double(*distance_metric)(T&, T&)){
// --------------------------------------------------------------------asigment 2
  #if DEBUG
  cout<<"assigment2\n";
  #endif
  vector<T*> *clusters=new vector<T*>[k];
  hash<T*> hasher;

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

  unordered_map<unsigned int, T*> left_to_classify;
  for(T& i : *data)//FIXME perini ora. mpori ta data na prepi na ine se set
    left_to_classify.insert(make_pair(hasher(&i),&i));

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

//TODO make update1 brute force
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
      if(tmp.second==DBL_MAX){//TODO what to do if lsh cant find options: brute force
        cerr<<"wrong lsh arguments!! no NN found!!\n\n";
        exit(1);
      }
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

my_curve get_curve_mean(unsigned int dimentions, vector<my_curve*> &cluster){//DEN prepi to cluster na ine 0!!
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

  my_curve *c2;
  vector<my_vector*>* a=new vector<my_vector*>[l];
  list<pair<unsigned int,unsigned int>>* ipairs;
  while(1){
    c2=c;
    for(unsigned int i=0;i<l;i++)
      a[i].clear();
    for(my_curve* it: cluster){
      ipairs=MinMatching(*c,*it);
      for(pair<unsigned int,unsigned int> ij: *ipairs)
        a[ij.first-1].push_back(it->vectors[ij.second-1]);
    }
    for(unsigned int i=0;i<l;i++)
      *c->vectors[i]=get_vector_mean(dimentions,a[i]);
    if(*c2==*c)
      break;
    //delete c2;
  }

  //delete c2;
  for(unsigned int i=0;i<l;i++)
    a[i].clear();
  delete[] a;

  return *c;
}

my_vector get_vector_mean(unsigned int dimentions, vector<my_vector*> &cluster){//DEN prepi to cluster na ine 0!!
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

double list_diff(vector<my_vector> *list1, vector<my_vector> *list2){
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

double list_diff(vector<my_curve> *list1, vector<my_curve> *list2){
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





template<class T>
list<double> * silhouette(vector<T*>* clusters,vector<T>* centers,unsigned int k,unsigned int n,double(*distance_metric)(T&, T&)){
  double *a= new double[n];
  double *b= new double[n];
  unsigned int clustersize = 0;
  unsigned int i = 0;//this is for a[i]
  unsigned int r = 0;//this is for b[r]
  double *clustercoef = new double[k];
  double **distarray;
  for (unsigned int x  = 0; x < n; x++) {
    b[x] = 0.0;
    a[x] = 0.0;
  }
  for (unsigned int x  = 0; x < k; x++) {
    clustercoef[x] = 0.0;
  }

  //NearestC & 2ndNearestC
  for(unsigned int j=0; j<k; j++){//take a center and its cluster(they sould be sync)

    clustersize = clusters[j].size();
    #if DEBUG
      std::cout << "j="<<j<<"k="<<k <<"i,r = "<<i<<","<<r <<"size="<<clustersize<<"with center"<<(*centers)[j].id<< '\n';
    #endif

    //create lower triangular matrix
    if (clustersize>1) {
      distarray = new double* [clustersize-1];//create distance array:for cluster
    }
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
    for (size_t l = 0; l < clustersize-1; l++) {
      delete distarray[l];
    }
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
    for (; i < sum+clustersize; i++) {
      // this is s(i)
      if(clustersize != 1){
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
      }else
        clustercoef[j]=0;
    }
    sum+=clustersize;
    clustercoef[j]=clustercoef[j]/(double)(clustersize);
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
  for (unsigned int i = 0; i < k; i++) {
    result->push_back(clustercoef[i]);
  }

  // delete a;
  // delete b;
  // delete clustercoef;

  return result;
}

template<class T>
list<double> *vectorizeandsilhouette(list<T*>* lclusters,list<T>* lcenters,unsigned int k,unsigned int n,double(*distance_metric)(T&, T&)){
  vector<T> centers(lcenters->begin(), lcenters->end());
  vector<T*> clusters;
  list<double> *result;
  for (T* c: (*lclusters)) {//naive approach
	   clusters.push_back(c);
	}

  // calls silhouette
  result = silhouette(&clusters,&centers,k,n,distance_metric);
  // clusters.clear();
  // centers.clear();

  return result;
}
