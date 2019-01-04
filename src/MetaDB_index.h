// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Nov 29, 2017
// Updated by Xiaoquan Su

#include <iostream>
#include <fstream>
#include <sstream>
#include "db.h"

#ifndef METADB_INDEX_H
#define METADB_INDEX_H

#define DEFAULT_INDEX_LEVEL 3

using namespace std;

class _MetaDB_Index{
      
      public:               
             _MetaDB_Index(){                                           
                              Set_Index_Level(DEFAULT_INDEX_LEVEL);                                                                                                                
                             }
             
             //For dev only
             _MetaDB_Index(int level){
                              Set_Index_Level(level);                                                                        
                             }
                                                                 
             int Get_IndexN();                                                 
             int Make_Index(float * abd_otu, float * abd_index);                        
             float Calc_Index(float * abd_1, float * abd_2, float t);   
             float Calc_Index_Cos(float * abd_1, float * abd_2, float t);   
             string Get_Index_Entry(int n);          
             
      private:              
             map <string, int> Index_hash;
             _PMDB Database;
             string Index_file;
             string Index_taxa_file;
             
             vector <string> Index_entry;
             vector <int> Order; //otu to index
             
             int IndexN;
             int OrderN;
             
             void Set_Index_Level(int level);
             void Init();             
             int Load_Index_Table();
             int Load_Index_Order();
            
      };

void _MetaDB_Index::Init(){                                    
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
    for (int i = 0; i < OrderN; i ++)
        if (Order[i] >= 0)
           abd_index[Order[i]] += abd_otu[i];
    
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

void _MetaDB_Index::Set_Index_Level(int level){ //For dev only
     
     switch(level){
               
               case 0: Index_file = Database.Get_Path() + "/index/Index_phylum.txt";
                       Index_taxa_file = Database.Get_Path() + "/index/Index_taxa_phylum.txt"; 
                       break;
               case 1: Index_file = Database.Get_Path() + "/index/Index_class.txt";
                       Index_taxa_file = Database.Get_Path() + "/index/Index_taxa_class.txt"; 
                       break;
               case 2: Index_file = Database.Get_Path() + "/index/Index_order.txt";
                       Index_taxa_file = Database.Get_Path() + "/index/Index_taxa_order.txt"; 
                       break;
               case 3: Index_file = Database.Get_Path() + "/index/Index_family.txt";
                       Index_taxa_file = Database.Get_Path() + "/index/Index_taxa_family.txt"; 
                       break;
               case 4: Index_file = Database.Get_Path() + "/index/Index_genus.txt";
                       Index_taxa_file = Database.Get_Path() + "/index/Index_taxa_genus.txt";
                       break;
               case 5: Index_file = Database.Get_Path() + "/index/Index_C73.txt";
                       Index_taxa_file = Database.Get_Path() + "/index/Index_C73_order.txt";
                       break;
               case 6: Index_file = Database.Get_Path() + "/index/Index_C75.txt";
                       Index_taxa_file = Database.Get_Path() + "/index/Index_C75_order.txt";
                       break;
               case 7: Index_file = Database.Get_Path() + "/index/Index_73.txt";
                       Index_taxa_file = Database.Get_Path() + "/index/Index_73_order.txt";
                       break;
               case 8: Index_file = Database.Get_Path() + "/index/Index_76.txt";
                       Index_taxa_file = Database.Get_Path() + "/index/Index_76_order.txt";
                       break;
               case 9: Index_file = Database.Get_Path() + "/index/Index_79.txt";
                       Index_taxa_file = Database.Get_Path() + "/index/Index_79_order.txt";
                       break;                             
              default: Index_file = Database.Get_Path() + "/index/Index_family.txt";
                       Index_taxa_file = Database.Get_Path() + "/index/Index_taxa_family.txt";
                       break;
               }
     Init();
     }

string _MetaDB_Index::Get_Index_Entry(int n){
       
       if (n < 0) return "NA";
       if (n >= IndexN) return "NA";
       
       return Index_entry[n];
       
       }
#endif
