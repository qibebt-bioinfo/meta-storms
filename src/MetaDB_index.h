// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Dec 20, 2019
// Updated by Xiaoquan Su

#include <iostream>
#include <fstream>
#include <sstream>
#include "db.h"
#include "math.h"
#ifndef METADB_INDEX_H
#define METADB_INDEX_H

#define DEFAULT_INDEX_LEVEL 3

using namespace std;

class _MetaDB_Index{
      
      public:               
             _MetaDB_Index(){}
             
             _MetaDB_Index(char db){
                              _PMDB database(db);                              
                              if (db == 'F') Database_path = database.Get_Func_Path();
                              else  Database_path = database.Get_Path();
                              Init();                                                                                                                                           
                             }
                                            
             int Get_IndexN();                                                 
             int Make_Index(float * abd_otu, float * abd_index);                        
             float Calc_Index(float * abd_1, float * abd_2, float t);   
             float Calc_Index_Cos(float * abd_1, float * abd_2, float t);   
             string Get_Index_Entry(int n);          
             
      private:              
             map <string, int> Index_hash;
             
             string Database_path;
             string Index_file;
             string Index_taxa_file;
             
             vector <string> Index_entry;
             vector <int> Order; //otu to index
             
             int IndexN;
             int OrderN;
             
             void Init();             
             int Load_Index_Table();
             int Load_Index_Order();
            
      };

void _MetaDB_Index::Init(){                                    
                  Index_file = Database_path + "/index/index.txt";
                  Index_taxa_file = Database_path + "/index/index_order.txt"; 
                  
                  IndexN = Load_Index_Table();
                  OrderN = Load_Index_Order();
                  }

int _MetaDB_Index::Load_Index_Table(){                                       
                    ifstream infile(Index_file.c_str(), ios::in);
                    if (!infile){
                                 cerr << "Error: Cannot open index file: " << Index_file << endl;
                                 return 0;
                                 }
                    string buffer;
                    while(getline(infile, buffer)){                                          
                                          if (buffer.size() == 0) continue;
                                          Index_entry.push_back(buffer);
                                          }
                    
                    infile.close();
                    infile.clear();
                                        
                    return Index_entry.size();
                    }

int _MetaDB_Index::Load_Index_Order(){                                       
                    ifstream infile(Index_taxa_file.c_str(), ios::in);
                    if (!infile){
                                 cerr << "Error: Cannot open index file: " << Index_taxa_file << endl;
                                 return 0;
                                 }
                    string buffer;
                    while(getline(infile, buffer)){                                          
                                          if (buffer.size() == 0) continue;
                                          int order;
                                          stringstream strin(buffer);
                                          strin >> order;
                                          Order.push_back(order);
                                          }
                    
                    infile.close();
                    infile.clear();
                    
                    return Order.size();
                    }

int _MetaDB_Index::Get_IndexN(){
    
    return IndexN;
    }

int _MetaDB_Index::Make_Index(float * abd_otu, float * abd_index){
    
    for (int i = 0; i < IndexN; i ++)
        abd_index[i] = 0;
    
    float index_sum = 0;
    
    for (int i = 0; i < OrderN; i ++)
        if (Order[i] >= 0){
           abd_index[Order[i]] += abd_otu[i];
           index_sum += abd_otu[i];
           }
    
    //norm
    index_sum /= 100.0;
    for (int i = 0; i < IndexN; i ++){
        abd_index[i] /= index_sum;
        }
    
    //count no zero
    int n_zero = 0;
    for (int i = 0; i < IndexN; i ++)
        if (abd_index[i] > 0) n_zero ++;
        
    return n_zero;
    }
 
float _MetaDB_Index::Calc_Index(float * abd_1, float * abd_2, float t){ // eu, smaller, better
      
      float index = 0;
      
      for (int i = 0; i < IndexN; i ++){
          
          float abd1 = (abd_1[i] > t) ? abd_1[i] : 0;
          float abd2 = (abd_2[i] > t) ? abd_2[i] : 0;
          
          index += pow(abd1 - abd2, 2);
          }                    
      return sqrt(index);
      }

/*
float _MetaDB_Index::Calc_Index(float * abd_1, float * abd_2, float t){ 
      
      return Calc_Index_Cos(abd_1, abd_2, t );
      
      }
 */
float _MetaDB_Index::Calc_Index_Cos(float * abd_1, float * abd_2, float t){ // cos, smaller, better
      
      float index = 0;
      
      float f_m_sum = 0;
      float f_n_sum = 0;
      float root_f = 0;
            
      for (int i = 0; i < IndexN; i ++){
          
          float abd1 = (abd_1[i] > t) ? abd_1[i] : 0;
          float abd2 = (abd_2[i] > t) ? abd_2[i] : 0;
          
          index += abd1 * abd2;
          
          f_m_sum += abd1 * abd1;
          f_n_sum += abd2 * abd2;
                   
          }
      
      
      root_f = (sqrt(f_m_sum * f_n_sum));
      
      if (root_f == 0) return 1;
      
      return index / root_f;
      }

string _MetaDB_Index::Get_Index_Entry(int n){
       
       if (n < 0) return "NA";
       if (n >= IndexN) return "NA";
       
       return Index_entry[n];
       
       }
#endif
