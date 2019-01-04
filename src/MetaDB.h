// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at June 14, 2018
// Updated by Xiaoquan Su
// comp.h
// Compatible with PM3 or above (16S)
// Incompatible for PM1, PM2 and 18S
// Added unweight
// Added HDD model

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>
#include <set>
#include <omp.h>
#include <stdio.h>

#include "comp.h"
#include "utility.h"
#include "MetaDB_index.h"
#include "table_format.h"
#include "qsort.h"
#include "hash.h"

#ifndef METADB_H
#define METADB_H

#define BUFF 10000000 
//string Database;
//string Coren/

class _MetaDB_Entry{
      public:
             friend class _MetaDB;
             
             _MetaDB_Entry(){
                            Abd = NULL;
                            Index_abd = NULL;
                            OTU_count = 0;
                            }
              
             _MetaDB_Entry(string id, float * abd, int abd_size, float * index_abd, int index_abd_size, string file){
                                  ID = id;
                                  File = file;
                                  OTU_count = 0;
                                  //abd
                                  if (abd_size > 0){
                                     Abd = new float [abd_size];
                                     memcpy(Abd, abd, abd_size * sizeof(float));
                                     }
                                  
                                  //index
                                  if (index_abd_size > 0){
                                     Index_abd = new float [index_abd_size];
                                     memcpy(Index_abd, index_abd, index_abd_size * sizeof(float));      
                                     }                    
                                  }
             /*
             ~_MetaDB_Entry(){
                              delete [] Abd;
                              delete [] Index_abd;
                              }
             */
             
             int Set_OTU_Count(float * abd, int abd_size){
                 //otu count
                 for (int i = 0; i < abd_size; i ++)
                                      OTU_count += (abd[i] > 0) ? 1: 0;
                 return OTU_count;
                 }
             
             int Get_OTU_Count(){
                 return OTU_count;
                 }
             
             string Get_ID(){
                    return ID;
                    }
             
             float * Get_Abd(){
                   return Abd;
                   }
             
             float * Get_Index_Abd(){
                   return Index_abd;
                   }
             string Get_File(){
                   return File;
                   }
             
      private:
              string ID;
              float * Abd;
              float * Index_abd;
              string File;   
              int OTU_count;                        
      };

class _MetaDB{
      
      public:
             _MetaDB(){
                       Is_cp_correct = false;
                       Mem_mode = 0;
                       Index_level = DEFAULT_INDEX_LEVEL;
                       };
                                                                          
             int Get_Sample_Number(){                 
                 return Samples.size();                 
                 }
                                            
             //Make database
             int Make_Database(const char * databasefile);        
                            
             //Load database
             int Load_Database(const char * databasefile, int coren);
             int Load_Database(const char * databasefile);
             int Load_Database_Level(const char * databasefile, int coren, int level); //for dev only
             int Load_Database_Level(const char * databasefile, int level); //for dev only
             
             //int Load_Database_Func(const char * databasefile);
             
             //Set prefix for HDD mode
             void Set_Sample_Hdd_Prefix(string p){
             	  Mem_mode = 1;
                  Sample_prefix = p;
                  }
                          
      protected:             
              vector <_MetaDB_Entry> Samples;  
              bool Is_cp_correct;
              int Index_level; //3: family
              int Mem_mode; //0: RAM; 1: HDD  
              string Sample_prefix; //prefix for HDD mode
                                          
             _Comp_Tree Comp_tree;
             _MetaDB_Index MetaDB_Index;   
             
             //_Comp_Tree_Func Comp_tree_func;     
			 
			 string Get_Sample_Hdd_Path(string id);                                    
      };

//Make database
int _MetaDB::Make_Database(const char * databasefile){
    
                            ofstream outfile(databasefile, ios::out);
                            if (!outfile){
                               cerr << "Error: Cannot open database file: " << databasefile << endl;
                               return 0;
                               }
                               
                           //output config
                            if (Is_cp_correct) outfile << "T";
                            else outfile << "F";
                            outfile << Index_level << endl; 
                           
                            for (int i = 0; i < Samples.size(); i ++){                                                                              
                                       
                                       //ID
                                       outfile << Samples[i].Get_ID() << endl;
                                       //abd
                                       outfile << "Abd:";
                                       for (int j = 0; j < Comp_tree.Get_LeafN(); j ++)
                                           if (Samples[i].Get_Abd()[j] > 0)
                                              outfile << "\t" << j << "\t" << Samples[i].Get_Abd()[j];
                                       outfile << endl;
                                       //index
                                       outfile << "Index:";
                                       for (int j = 0; j < MetaDB_Index.Get_IndexN(); j ++)
                                           if (Samples[i].Get_Index_Abd()[j] > 0)
                                              outfile << "\t" << j << "\t" << Samples[i].Get_Index_Abd()[j];
                                       outfile << endl;
                                       //Other
                                       outfile << Samples[i].Get_File() << endl;
                                                                
                                }
                            
                            return Samples.size();
    }

//Load
/*
int _MetaDB::Load_Database(const char * databasefile, int coren){				  	
                 FILE * fptr = fopen(databasefile, "r");
                 
                 if (fptr == NULL){
                              cerr << "Error: Cannot open database file: " << databasefile << endl;
                              return 0;
                              }
                                                   
                 char * str_buffer = new char [BUFF];
                 vector <string> file_buffer;
                 while(fgets(str_buffer, BUFF-1, fptr) != NULL){
                                        
                                        if (str_buffer[strlen(str_buffer)-1] == '\n')
												str_buffer[strlen(str_buffer)-1] = '\0';  
                                        string a_line = str_buffer;
                                        file_buffer.push_back(a_line);
                                        
                                        }
                 fclose(fptr); 
                 delete [] str_buffer;
				 //debug
                 cout << "IO finished" << endl;
                 
                 int iter = 0;   
                 
                 //config                                               
                 while(iter < file_buffer.size()){ //load config
                                       string buffer = file_buffer[iter];
                                       if (buffer.size() == 0) continue;
                                       if (buffer[0] == '#') continue;
                                       if (buffer[0] == 'T')
                                                     Is_cp_correct = true;
                                       else Is_cp_correct = false;
                                       iter ++;              
                                       break;
                                       }
                 
                 long sample_size = (file_buffer.size() - iter) / 4;
                 //debug
                 cout << "Samples: " << sample_size << endl;
                 
                 Samples = vector <_MetaDB_Entry> (sample_size, _MetaDB_Entry());
                 
                 omp_set_num_threads(coren);
                 #pragma omp parallel for schedule(dynamic, 1)                                                            
                 for (int i = 0; i < sample_size; i ++){ //load samples
                                                                              
                                       string id = file_buffer[i * 4 + iter];
                                       string abd_buffer = file_buffer[i * 4 + 1 + iter];
                                       string index_abd_buffer = file_buffer[i * 4 + 2 + iter];
                                       string file = file_buffer[i * 4 + 3 + iter];
                                       string temp;
                                       
                                       //abd                                       
                                       float * abd = new float [Comp_tree.Get_LeafN()];
                                       memset(abd, 0, Comp_tree.Get_LeafN() * sizeof(float));
                                       if (Mem_mode ==0){ //RAM mode                                                                                 
                                          stringstream strin_abd(abd_buffer);
                                       
                                          int a_abd_pos;
                                          float a_abd;
                                       
                                          strin_abd >> temp;                                                                            
                                          while(strin_abd >> a_abd_pos){
                                                       strin_abd >> a_abd;
                                                       abd[a_abd_pos] = a_abd;
                                                       }
                                       }
                                       //index_abd
                                       float * index_abd = new float [MetaDB_Index.Get_IndexN()];
                                       memset(index_abd, 0, MetaDB_Index.Get_IndexN() * sizeof(float));
                                       
                                       MetaDB_Index.Make_Index(abd, index_abd);
                                                                              
                                       //add sample
                                       if (Mem_mode == 0){ //RAM mode
                                          Samples[i] = _MetaDB_Entry(id, abd, Comp_tree.Get_LeafN(), index_abd, MetaDB_Index.Get_IndexN(), file);
                                          Samples[i].Set_OTU_Count(abd, Comp_tree.Get_LeafN());     
                                      	  }
                                       else //HDD mode
                                            Samples[i] = _MetaDB_Entry(id, NULL, 0, index_abd, MetaDB_Index.Get_IndexN(), file);                                                                              
                                       
                                       delete [] abd;
                                       delete [] index_abd;
                                       }                                  
                 file_buffer.clear();
                 vector <string> ().swap(file_buffer);
                 
                 return Samples.size();
                 }
*/
int _MetaDB::Load_Database(const char * databasefile, int coren){
				
                 FILE * fptr = fopen(databasefile, "r");                 
                 if (fptr == NULL){
                              cerr << "Error: Cannot open database file: " << databasefile << endl;
                              return 0;
                              }
                 //config                                               
                //while(iter < file_buffer.size()){ //load config
                char * file_buffer = new char [BUFF];
                while(fgets(file_buffer, BUFF-1, fptr) != NULL){
                                       string buffer = file_buffer;
                                       if (buffer.size() == 0) continue;
                                       if (buffer[0] == '#') continue;
                                       if (buffer[0] == 'T')
                                                     Is_cp_correct = true;
                                       else Is_cp_correct = false;            
                                       break;
                                       }               
         
				while(fgets(file_buffer, BUFF-1, fptr) != NULL){
                                       if (file_buffer[strlen(file_buffer)-1] == '\n')
												file_buffer[strlen(file_buffer)-1] = '\0';                                       
                                       string id = file_buffer;
                                       fgets(file_buffer, BUFF-1, fptr);
                                       string abd_buffer = file_buffer;
                                       fgets(file_buffer, BUFF-1, fptr);
                                       string index_abd_buffer = file_buffer;
                                       fgets(file_buffer, BUFF-1, fptr);
                                       string file = "";
                                       string temp;
                                       
                                       //abd                                       
                                       float * abd = new float [Comp_tree.Get_LeafN()];
                                       memset(abd, 0, Comp_tree.Get_LeafN() * sizeof(float));
                                       if (Mem_mode ==0){ //RAM mode                                                                                 
                                          stringstream strin_abd(abd_buffer);
                                       
                                          int a_abd_pos;
                                          float a_abd;
                                       
                                          strin_abd >> temp;                                                                            
                                          while(strin_abd >> a_abd_pos){
                                                       strin_abd >> a_abd;
                                                       abd[a_abd_pos] = a_abd;
                                                       } 
                                       }
                                       //index_abd
                                       float * index_abd = new float [MetaDB_Index.Get_IndexN()];
                                       memset(index_abd, 0, MetaDB_Index.Get_IndexN() * sizeof(float));
                                       stringstream strin_index(index_abd_buffer);
                                                                              
                                       int a_index_pos;
                         
						               float a_index;
                                       
                                       strin_index >> temp;
                                       while(strin_index >> a_index_pos){
                                                         strin_index >> a_index;
                                                         index_abd[a_index_pos] = a_index;
                                                         }
									                                         
                                       //add sample
                                       _MetaDB_Entry a_sample;
                                       if (Mem_mode == 0){ //RAM mode
                                          a_sample = _MetaDB_Entry(id, abd, Comp_tree.Get_LeafN(), index_abd, MetaDB_Index.Get_IndexN(), file);
                                          a_sample.Set_OTU_Count(abd, Comp_tree.Get_LeafN());  
                                      	}
                                       else //HDD mode
                                            a_sample = _MetaDB_Entry(id, NULL, 0, index_abd, MetaDB_Index.Get_IndexN(), file);
                                       
                                       Samples.push_back(a_sample);
									                                                                                    
                                       delete [] abd;
                                       delete [] index_abd;  
									                                                                              
                                }
				 fclose(fptr);
				 
                 return Samples.size();
                 }

int _MetaDB::Load_Database(const char * databasefile){
		
		Load_Database(databasefile, 1);
		
		}
/*
int _MetaDB::Load_Database_Level(const char * databasefile, int coren, int level){ //for dev only
				Index_level = level;
				  	
                 FILE * fptr = fopen(databasefile, "r");
                 
                 if (fptr == NULL){
                              cerr << "Error: Cannot open database file: " << databasefile << endl;
                              return 0;
                              }
                                                   
                 char * str_buffer = new char [BUFF];
                 vector <string> file_buffer;
                 while(fgets(str_buffer, BUFF-1, fptr) != NULL){
                                        
                                        if (str_buffer[strlen(str_buffer)-1] == '\n')
												str_buffer[strlen(str_buffer)-1] = '\0';  
                                        string a_line = str_buffer;
                                        file_buffer.push_back(a_line);
                                        
                                        }
                 fclose(fptr); 
                 delete [] str_buffer;
				 //debug
                 cout << "IO finished" << endl;
                 
                 int iter = 0;   
                 
                 //config                                               
                 while(iter < file_buffer.size()){ //load config
                                       string buffer = file_buffer[iter];
                                       if (buffer.size() == 0) continue;
                                       if (buffer[0] == '#') continue;
                                       if (buffer[0] == 'T')
                                                     Is_cp_correct = true;
                                       else Is_cp_correct = false;
                                       iter ++;              
                                       break;
                                       }
                 
                
				 MetaDB_Index = _MetaDB_Index(Index_level);
                 
                 long sample_size = (file_buffer.size() - iter) / 4;
                 //debug
                 cout << "Samples: " << sample_size << endl;
                 
                 Samples = vector <_MetaDB_Entry> (sample_size, _MetaDB_Entry());
                 
                 omp_set_num_threads(coren);
                 #pragma omp parallel for schedule(dynamic, 1)                                                            
                 for (int i = 0; i < sample_size; i ++){ //load samples
                                                                              
                                       string id = file_buffer[i * 4 + iter];
                                       string abd_buffer = file_buffer[i * 4 + 1 + iter];
                                       string index_abd_buffer = file_buffer[i * 4 + 2 + iter];
                                       string file = file_buffer[i * 4 + 3 + iter];
                                       string temp;
                                       
                                       //abd                                       
                                       float * abd = new float [Comp_tree.Get_LeafN()];
                                       memset(abd, 0, Comp_tree.Get_LeafN() * sizeof(float));
                                       if (Mem_mode ==0){ //RAM mode                                                                                 
                                          stringstream strin_abd(abd_buffer);
                                       
                                          int a_abd_pos;
                                          float a_abd;
                                       
                                          strin_abd >> temp;                                                                            
                                          while(strin_abd >> a_abd_pos){
                                                       strin_abd >> a_abd;
                                                       abd[a_abd_pos] = a_abd;
                                                       }
                                       }
                                       //index_abd
                                       float * index_abd = new float [MetaDB_Index.Get_IndexN()];
                                       memset(index_abd, 0, MetaDB_Index.Get_IndexN() * sizeof(float));
                                       
                                       MetaDB_Index.Make_Index(abd, index_abd);
                                                                              
                                       //add sample
                                       if (Mem_mode == 0){ //RAM mode
                                          Samples[i] = _MetaDB_Entry(id, abd, Comp_tree.Get_LeafN(), index_abd, MetaDB_Index.Get_IndexN(), file);
                                          Samples[i].Set_OTU_Count(abd, Comp_tree.Get_LeafN());     
                                      	  }
                                       else //HDD mode
                                            Samples[i] = _MetaDB_Entry(id, NULL, 0, index_abd, MetaDB_Index.Get_IndexN(), file);                                                                              
                                       
                                       delete [] abd;
                                       delete [] index_abd;
                                       }                                  
                 file_buffer.clear();
                 vector <string> ().swap(file_buffer);
                 
                 return Samples.size();
                 }
*/
int _MetaDB::Load_Database_Level(const char * databasefile, int coren, int level){//for dev only 
				 Index_level = level;
				 MetaDB_Index = _MetaDB_Index(Index_level);
				 
                 FILE * fptr = fopen(databasefile, "r");                 
                 if (fptr == NULL){
                              cerr << "Error: Cannot open database file: " << databasefile << endl;
                              return 0;
                              }
                 //config                                               
                //while(iter < file_buffer.size()){ //load config
                char * file_buffer = new char [BUFF];
                while(fgets(file_buffer, BUFF-1, fptr) != NULL){
                                       string buffer = file_buffer;
                                       if (buffer.size() == 0) continue;
                                       if (buffer[0] == '#') continue;
                                       if (buffer[0] == 'T')
                                                     Is_cp_correct = true;
                                       else Is_cp_correct = false;            
                                       break;
                                       }               
         
				while(fgets(file_buffer, BUFF-1, fptr) != NULL){
                                       if (file_buffer[strlen(file_buffer)-1] == '\n')
												file_buffer[strlen(file_buffer)-1] = '\0';                                          
                                       string id = file_buffer;
                                       fgets(file_buffer, BUFF-1, fptr);
                                       string abd_buffer = file_buffer;
                                       fgets(file_buffer, BUFF-1, fptr);
                                       string index_abd_buffer = file_buffer;
                                       fgets(file_buffer, BUFF-1, fptr);
                                       string file = "";
                                       string temp;
                                       
                                       //abd                                       
                                       float * abd = new float [Comp_tree.Get_LeafN()];
                                       memset(abd, 0, Comp_tree.Get_LeafN() * sizeof(float));
                                       if (Mem_mode ==0){ //RAM mode                                                                                 
                                          stringstream strin_abd(abd_buffer);
                                       
                                          int a_abd_pos;
                                          float a_abd;
                                       
                                          strin_abd >> temp;                                                                            
                                          while(strin_abd >> a_abd_pos){
                                                       strin_abd >> a_abd;
                                                       abd[a_abd_pos] = a_abd;
                                                       } 
                                       }
                                       //index_abd
                                       float * index_abd = new float [MetaDB_Index.Get_IndexN()];
                                       memset(index_abd, 0, MetaDB_Index.Get_IndexN() * sizeof(float));
                                       
                                       MetaDB_Index.Make_Index(abd, index_abd);
									                                         
                                       //add sample
                                       _MetaDB_Entry a_sample;
                                       if (Mem_mode == 0){ //RAM mode
                                          a_sample = _MetaDB_Entry(id, abd, Comp_tree.Get_LeafN(), index_abd, MetaDB_Index.Get_IndexN(), file);
                                          a_sample.Set_OTU_Count(abd, Comp_tree.Get_LeafN());  
                                      	}
                                       else //HDD mode
                                            a_sample = _MetaDB_Entry(id, NULL, 0, index_abd, MetaDB_Index.Get_IndexN(), file);
                                       
                                       Samples.push_back(a_sample);
									                                                                                    
                                       delete [] abd;
                                       delete [] index_abd;  
									                                                                              
                                }
				 fclose(fptr);
				 
                 return Samples.size();
                 }
int _MetaDB::Load_Database_Level(const char * databasefile, int level){ //for dev only
		
		Load_Database_Level(databasefile, 1, level);
		
		}

string _MetaDB::Get_Sample_Hdd_Path(string id){
	
	return (Sample_prefix + "/" + id + ".dat");
	
} 
                          
//Load_func 
/*
int _MetaDB::Load_Database_Func(const char * databasefile){
    
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
                                       //if (Mem_mode == 0) //RAM mode
                                          a_sample = _MetaDB_Entry(id, abd, Comp_tree_func.Get_GeneN(), NULL, 0, file);
                                       //else //HDD mode
                                       //     a_sample = _MetaDB_Entry(id, NULL, 0, NULL, 0, file);
                                       
                                       a_sample.Set_OTU_Count(abd, Comp_tree_func.Get_GeneN());      
                                       
                                       Samples.push_back(a_sample); 
                                        
                                       delete [] abd;
                                       //delete [] index_abd;
                                       }
                 
                 infile.close();
                 infile.clear();
                 
                 return Samples.size();
                 }
*/
   
#endif
