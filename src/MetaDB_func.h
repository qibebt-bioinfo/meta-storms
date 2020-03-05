// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Dec 11, 2018
// Updated by Xiaoquan Su
// MetaDB__func

#include <iostream>
#include "MetaDB.h"
#include "comp_sam_func.h"

#ifndef METADB_FUNC_H
#define METADB_FUNC_H

class _MetaDB_func : public _MetaDB{
      
      public:
             int Load_Database_Func(const char * databasefile);    
                          
      protected:               
             _Comp_Tree_Func Comp_tree_func;     
                                  
      };


//Load_func 

int _MetaDB_func::Load_Database_Func(const char * databasefile){
    
                 ifstream infile(databasefile, ios::in);
                 if (!infile){
                              cerr << "Error: Cannot open database file: " << databasefile << endl;
                              return 0;
                              }
                 
                 string buffer;
                 while(getline(infile, buffer)){ //load config
                                       if (buffer.size() == 0) continue;
                                       if (buffer[0] == '#') continue;
                                       if (buffer[0] == 'T')
                                                     Is_cp_correct = true;
                                       else Is_cp_correct = false;
                                                     
                                       if (buffer.size() > 1)
                                          Index_level = buffer[1] - '0';
                                       break;
                                       }
                                                                              
                 while(getline(infile, buffer)){//load samples
                                       
                                       //cout << "Loading sample: " << Samples.size() + 1 << endl; //debug
                                       
                                       if (buffer.size() == 0) continue;
                                       if (buffer[0] == '#') continue;
                                       string id;
                                       string file;
                                       string abd_buffer;
                                       string index_abd_buffer;
                                       string temp;
                                       id = buffer;
                                       //abd
                                       getline(infile, abd_buffer);                                       
                                       float * abd = new float [Comp_tree_func.Get_GeneN()];
                                       memset(abd, 0, Comp_tree_func.Get_GeneN() * sizeof(float));
                                       stringstream strin_abd(abd_buffer);
                                       
                                       int a_abd_pos;
                                       float a_abd;
                                       
                                       strin_abd >> temp;                                                                            
                                       while(strin_abd >> a_abd_pos){
                                                       strin_abd >> a_abd;
                                                       abd[a_abd_pos] = a_abd;
                                                       }
                     
                                       //index_abd
                                       //no index for func
                                       getline(infile, index_abd_buffer);
                                       
                                       //file path
                                       getline(infile, file);
                                       
                                       //add sample
                                       _MetaDB_Entry a_sample;
                                       
                                       //add sample
                                       if (Mem_mode == 0){ //RAM mode
                                          a_sample = _MetaDB_Entry(id, abd, Comp_tree_func.Get_GeneN(), NULL, 0, file);
                                          a_sample.Set_OTU_Count(abd, Comp_tree_func.Get_GeneN());
                                          }
                                       else //HDD mode
                                            a_sample = _MetaDB_Entry(id, NULL, 0, NULL, 0, file);
                                                                                                                           
                                       Samples.push_back(a_sample); 
                                        
                                       delete [] abd;
                                       //delete [] index_abd;
                                       }
                 
                 infile.close();
                 infile.clear();
                 
                 return Samples.size();
                 }
   
#endif
