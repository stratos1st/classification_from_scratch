#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <stdlib.h>

#include <unistd.h>
#include <string.h>

#include "util.hpp"
#include "my_vector.hpp"
#include "my_curve.hpp"
#include "lsh.hpp"
#include "clustering_funcs.hpp"

using namespace std;

#define DEBUG 0

template<class T>
void run_algorithms(unsigned int k, unsigned int max_iterations, unsigned int lsh_window, unsigned int g_no,
                    unsigned int grids_no, unsigned int container_sz, unsigned int lsh_l, unsigned int max_curve_sz, double grid_delta,
                    bool vector_input, bool stop_when_centers, bool brute_update_1, short int complete_flag, double center_tol, double pad_number,
                    char input_file[], char out_file[], char options_file[], list<T>* (*read_file)(string,unsigned int),
                    double (*distance_metric)(T&,T&));
template void run_algorithms<>(unsigned int k, unsigned int max_iterations, unsigned int lsh_window, unsigned int g_no,
                    unsigned int grids_no, unsigned int container_sz, unsigned int lsh_l, unsigned int max_curve_sz, double grid_delta,
                    bool vector_input, bool stop_when_centers, bool brute_update_1, short int complete_flag, double center_tol, double pad_number,
                    char input_file[], char out_file[], char options_file[], list<my_vector>* (*read_file)(string,unsigned int),
                    double (*distance_metric)(my_vector&,my_vector&));
template void run_algorithms<>(unsigned int k, unsigned int max_iterations, unsigned int lsh_window, unsigned int g_no,
                    unsigned int grids_no, unsigned int container_sz, unsigned int lsh_l, unsigned int max_curve_sz, double grid_delta,
                    bool vector_input, bool stop_when_centers, bool brute_update_1, short int complete_flag, double center_tol, double pad_number,
                    char input_file[], char out_file[], char options_file[], list<my_curve>* (*read_file)(string,unsigned int),
                    double (*distance_metric)(my_curve&,my_curve&));

int main(int argc, char** argv){
  srand (time(NULL));
  char input_file[100]("./Input/DataVectors_5_500x100.csv"),out_file[100]("my.out"),options_file[100]("cluster.conf");
  unsigned int k=4, max_iterations=50,lsh_window=6000,g_no=4,grids_no=4,container_sz=10,lsh_l=4,max_curve_sz=10;
  bool vector_input=true,stop_when_centers=true,brute_update_1=false;
  int complete_flag=0;// ine bool apla den iparxi atoi gia bool ke barieme na balo strtoul :'(
  double center_tol=0.01,pad_number=99999.99999,grid_delta=0.01;
  my_curve::curve_tol=0.001;
  my_vector::vector_tol=0.01;

  cout<<"Classification starting!!\n";

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
        complete_flag=atoi(optarg);
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
      sscanf(str,"brute_update_1: %u",&brute_update_1);
      sscanf(str,"lsh_curves_grids_delta: %u",&grid_delta);
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
  <<"\n\tmax_curve_size= "<<max_curve_sz<<"\n\tgrid_delta= "<<grid_delta<<"\n\tbrute_update_1= "<<brute_update_1
  <<"\n\tstop_when_centers_dont_change= "<<stop_when_centers<<"\n\tcurve_tolerance= "<<my_curve::curve_tol
  <<"\n\tvector_tolerance= "<<my_vector::vector_tol<<"\n\tcenter_tolerance= "<<center_tol
  <<"\n\tmax_iterations= "<<max_iterations<<"\n\tinput_file= "<<input_file
  <<"\n\tout_file= "<<out_file<<"\n\tcomplete_flag= "<<complete_flag<<endl;

  if(vector_input)
    run_algorithms(k, max_iterations, lsh_window, g_no,
                        grids_no, container_sz, lsh_l, max_curve_sz, grid_delta,
                        vector_input, stop_when_centers, brute_update_1, complete_flag, center_tol, pad_number,
                        input_file, out_file, options_file, read_vector_file, manhattan_distance);
  else
    run_algorithms(k, max_iterations, lsh_window, g_no,
                        grids_no, container_sz, lsh_l, max_curve_sz, grid_delta,
                        vector_input, stop_when_centers, brute_update_1, complete_flag, center_tol, pad_number,
                        input_file, out_file, options_file, read_curve_file, Dtw);

  cout<<"\nClassification ended!!\n";
  return 0;
}

template<class T>
void run_algorithms(unsigned int k, unsigned int max_iterations, unsigned int lsh_window, unsigned int g_no,
                    unsigned int grids_no, unsigned int container_sz, unsigned int lsh_l, unsigned int max_curve_sz, double grid_delta,
                    bool vector_input, bool stop_when_centers, bool brute_update_1, short int complete_flag, double center_tol, double pad_number,
                    char input_file[], char out_file[], char options_file[], list<T>* (*read_file)(string,unsigned int),
                    double (*distance_metric)(T&,T&)){
  short int function_matrix[3];

  ofstream fout(out_file);
  if (!fout.is_open()) {
    cerr<<"\n\n!! problem creating out file!!\n\n";
    exit(1);
  }

  //-----------------------------------------------------read input file
  list <T>* data_tmp=read_file(input_file,0);

  cout<<"reading done\n";
  #if DEBUG
  cout<<"data\n";
  for(auto i : *data_tmp)
    i.print_vec();
  #endif

  lsh *lsh_model;
  if(vector_input)
    lsh_model=new lsh_vector(data_tmp->front().get_dimentions(),lsh_l,lsh_window,g_no,container_sz);
  else{
    GridHash::delta = grid_delta;
    lsh_model=new lsh_curve(data_tmp->front().get_dimentions(),max_curve_sz,grids_no,lsh_window,g_no,pad_number,container_sz);
  }

  lsh_model->train(data_tmp);
  cout<<"lsh training done\n";

  vector<T>* data=new vector<T>;//to metatrepo se vector gt i palia sinartisi diabasmatos epestrefe list
  for(auto i: *data_tmp)
    data->push_back(i);
  data_tmp->clear();
  delete data_tmp;

  vector<T> *old_centers,*centers;
  vector<T*> *old_clusters,*clusters;
  unsigned int max_iterations_const=max_iterations;
  using namespace std::chrono;

  for(function_matrix[0]=0;function_matrix[0]<=1;function_matrix[0]++){
    for(function_matrix[1]=0;function_matrix[1]<=1;function_matrix[1]++){
      for(function_matrix[2]=0;function_matrix[2]<=1;function_matrix[2]++){
        //start clock
        auto start = high_resolution_clock::now();

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
          if(!brute_update_1)
            cout<<"update1\n";
          else
            cout<<"update1_brute\n";
        else
          cout<<"update2\n";
        cout<<endl;

        //-----------------------------------------------------------------------initialization
        if(function_matrix[0])
          centers=initialization1(data,k);
        else
          centers=initialization2(data,k,distance_metric);

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
          if(function_matrix[2]){
            if(!brute_update_1)
              centers=update1(data, centers, clusters,k,lsh_model,get_mean,distance_metric);
            else
              centers=update1_brute(data, centers, clusters,k,distance_metric);
          }
          else
            centers=update2(data, centers, clusters,k,get_mean);

          // sin8iki break:: nea kentra konta sta palia
          if(stop_when_centers)
            if(old_centers_equal_new_centers(old_centers,centers,center_tol,distance_metric))
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

          iteration++;
          //sin8iki break:: max iterations reached
          if(max_iterations==1){
            cout<<"max iterations reached\n";
            break;
          }
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

        //stop clock
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<seconds>(stop - start);

        //calculate silhouette
        list<double> *silhouette_vals = silhouette(clusters,centers,centers->size(),data->size(),distance_metric);

        ////--------------------------------------------print to out file
        if (fout.is_open()){
          fout<<"Algorithm: ";
          if(function_matrix[0])
            fout<<"initialization1 + ";
          else
            fout<<"initialization2 + ";
          if(function_matrix[1])
            fout<<"assigment1 + ";
          else
            fout<<"assigment2 + ";
          if(function_matrix[2])
            if(!brute_update_1)
              fout<<"update1\n";
            else
              fout<<"update1_brute\n";
          else
            fout<<"update2\n";
          for(unsigned int i=0;i<k;i++){
            fout<<"CLUSTER-"<<i+1<<"{ size: "<<clusters[i].size()<<", centroid: ";
            if(function_matrix[2])
              fout<<centers->at(i).id<<" }\n";
            else{
              streambuf *coutbuf = std::cout.rdbuf();
              cout.rdbuf(fout.rdbuf());//redirect cout to file
              centers->at(i).print_vec();
              cout.rdbuf(coutbuf);//revert
              fout<<" }\n";
            }
          }
          fout<<"clustering_time: "<<duration.count()<<endl;
          fout<<"Silhouette: [ ";
          for(auto i: *silhouette_vals)
            fout<<i<<", ";
          fout<<"]\n";
          if(complete_flag){
            for(unsigned int i=0;i<k;i++){
              fout<<"CLUSTER-"<<i+1<<"{ ";
              for(auto j: clusters[i])
                fout<<j->id<<",";
              fout<<" }\n";
            }
          }
          fout<<"\n\n"<<flush;
        }
        else{
          cerr<<"\n\n!! problem writing to out file!!\n\n";
          exit(1);
        }

        //--------------------------------------------deletes
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
        silhouette_vals->clear();
        delete silhouette_vals;

      }
    }
  }

  fout.close();
  delete lsh_model;
  data->clear();
  delete data;
}
