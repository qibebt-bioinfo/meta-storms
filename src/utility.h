// Updated at Dec 15, 2017
// Updated by Xiaoquan Su
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS

#ifndef _UTILITY_H
#define _UTILITY_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include <stdlib.h>
#include <string.h>
#include <sys/dir.h>
#include <sys/stat.h>
#include <unistd.h>

#define BUFFER_SIZE 5000

#include "hash.h"
using namespace std;

string Check_Env(){
    
    
    if (getenv("MetaStorms") == NULL){
                               
                               cerr << "Error: Please set the environment variable \"MetaStorms\" to the directory" << endl;
                               exit(0);
                               
                                   }
    
    string path =  getenv("MetaStorms");
    return path;
    }

int Check_Path(const char * path, int type){
    
    if (strlen(path) < 1) return 0;
    
    DIR *pDir = opendir(path);
            
    if(pDir!=NULL){
                  closedir(pDir);                   
                  if (type == 0){
                     string command = "rm -rf ";
                     command += path;
                     system(command.c_str());
                     mkdir(path, 0755);                          
                     }
                  }                  
    else 
         mkdir(path, 0755);           
    return 0;
    
    }

bool Check_Path(const char * path){
    
    if (strlen(path) < 1) return false;
    
    DIR *pDir = opendir(path);
    
    if (pDir != NULL){
             closedir(pDir);
             return true;
             }
    
    return false;
    }

bool Check_File(const char * file){
    
    fstream infile(file, ifstream::in);
    
    if (!infile){
                 
                 cerr << "Error: Cannot open file : " << file << endl;
                 return false;
                 }
    
    infile.close();
    infile.clear();
    
    return true;
    
    }

string Check_OTU(string otu){
       
       string a_otu = otu;
            if (a_otu.size() > 4 ){
                               string prefix_4 = a_otu.substr(0, 4);
                               if (( prefix_4 == "otu_") || ( prefix_4 == "OTU_"))
                                     a_otu = a_otu.substr(4, a_otu.size() - 4);
                               }
       return a_otu;
       }

unsigned int Get_Count(const char * infilename){
         
         ifstream infile(infilename, ifstream::in);
         
         if (!infile){
                      
                      cerr << "Error: Cannot open file : " << infilename << endl;
                      return 0;
                      
                      }
         
         string buffer;
         unsigned int count = 0;
         
         while (getline(infile, buffer)){
               
               if (buffer[0] == '>') count ++;
               
               }
         
         infile.close();
         infile.clear();
         
         return count;
         }

int Load_Copy_Number(string database_path, hash_map <string, float, std_string_hash> & cp_number){
    
    string copy_number_file = database_path + "16s_copy_number.txt";

    ifstream in_cp_number(copy_number_file.c_str(), ifstream::in);
    if (!in_cp_number){
                       cerr << "Error: Cannot open copy number table : " << copy_number_file << endl;
                       return 0;
                       } 
    
    int count = 0;
    string buffer;
    getline(in_cp_number, buffer);
    while(getline(in_cp_number, buffer)){
                                
                                stringstream strin(buffer);
                                string id;
                                float cp_no;
                                strin >> id >> cp_no;
                                
                                if (cp_number.count(id) == 0)
                                                           cp_number[id] = cp_no;
                                else cerr << "Warning, dup id : " << id << endl;
                                
                                count ++;
                                }
    
    in_cp_number.close();
    in_cp_number.clear();
    
    return count;
    }


int Load_ID(const char * idfilename, vector <string> & ID){
    
    ifstream infile(idfilename, ifstream::in);
             if (!infile){
                       cerr << "Error: Cannot open ID file: " << idfilename << endl;
                       return 0;
                       }
    
    string buffer;
    int count = 0;
    
    while(getline(infile, buffer)){
                          
                          if (buffer.size() == 0) continue;
                          
                          stringstream strin(buffer);
                          string id;
                          strin >> id;
                          
                          ID.push_back(id);
                          
                          count ++;
                          
                          }
    
    infile.close();
    infile.clear();        
    
    return count;
    }

int Load_List(const char * listfilename, vector <string> & list){
    
    ifstream infile(listfilename, ifstream::in);
             if (!infile){
                       cerr << "Error: Cannot open list file: " << listfilename << endl;
                       return 0;
                       }
    
    string buffer;
    int count = 0;
    
    while(getline(infile, buffer)){
                          
                          if (buffer.size() == 0) continue;
                          
                          string temp;
                          stringstream strin(buffer);
                          
                          while(strin >> temp);
                                                    
                          list.push_back(temp);
                          
                          count ++;
                          
                          }
    
    infile.close();
    infile.clear();        
    
    return count;
    }

int Load_List(const char * listfilename, vector <string> & list, string prefix){ //with prefix
    
    ifstream infile(listfilename, ifstream::in);
             if (!infile){
                       cerr << "Error: Cannot open list file: " << listfilename << endl;
                       return 0;
                       }
    
    string buffer;
    int count = 0;
    
    while(getline(infile, buffer)){
                          
                          if (buffer.size() == 0) continue;
                          
                          string temp;
                          stringstream strin(buffer);
                          
                          while(strin >> temp);
                                                    
                          list.push_back(prefix + temp);
                          
                          count ++;
                          
                          }
    
    infile.close();
    infile.clear();        
    
    return count;
    }

int Load_List(const char * listfilename, vector <string> & list, vector <string> & ids){//with id
    
    ifstream infile(listfilename, ifstream::in);
             if (!infile){
                       cerr << "Error: Cannot open list file: " << listfilename << endl;
                       return 0;
                       }
    
    string buffer;
    int count = 0;
    
    while(getline(infile, buffer)){
                          
                          if (buffer.size() == 0) continue;
                          
                          vector <string> str_temp;
                          string temp;
                          stringstream strin(buffer);
                          
                          while(strin >> temp)
                                      str_temp.push_back(temp);
                          
                          if (str_temp.size() >= 2){
                                              
                                              ids.push_back(str_temp[0]);
                                              list.push_back(str_temp[1]);
                                              
                                              }
                          
                          else{                               
                               list.push_back(buffer);   
                               int name_begin = buffer.find_first_of('/');
                               int name_end = buffer.find_last_of('/');
                               if (name_begin == name_end)
                                              name_begin = 0;
                               else name_begin = buffer.rfind('/', name_end-1)+1;
                               ids.push_back(buffer.substr(name_begin, name_end - name_begin));
                               }
                          count ++;
                          
                          }
    
    infile.close();
    infile.clear();        
    
    return count;
    }

int Load_List(const char * listfilename, vector <string> & list, vector <string> & ids, string prefix){//with id
    
    ifstream infile(listfilename, ifstream::in);
             if (!infile){
                       cerr << "Error: Cannot open list file: " << listfilename << endl;
                       return 0;
                       }
    
    string buffer;
    int count = 0;
    
    while(getline(infile, buffer)){
                          
                          if (buffer.size() == 0) continue;
                          
                          vector <string> str_temp;
                          string temp;
                          stringstream strin(buffer);
                          
                          while(strin >> temp)
                                      str_temp.push_back(temp);
                          
                          if (str_temp.size() >= 2){
                                              
                                              ids.push_back(str_temp[0]);
                                              list.push_back(prefix + str_temp[1]);
                                              
                                              }
                          
                          else{                               
                               list.push_back(prefix + buffer);   
                               int name_begin = buffer.find_first_of('/');
                               int name_end = buffer.find_last_of('/');
                               if (name_begin == name_end)
                                              name_begin = 0;
                               else name_begin = buffer.rfind('/', name_end-1)+1;
                               ids.push_back(buffer.substr(name_begin, name_end - name_begin));
                               }
                          count ++;
                          
                          }
    
    infile.close();
    infile.clear();        
    
    return count;
    }

void Get_Random(int n, int * order, int s, int loop){
    
    srand((int)time(NULL) + loop * 3000);
    
    int * order_table = new int [n];
    
    for(int i = 0; i< n; i++){
            order_table[i] = i;
            }
    for (int i = 0; i < s; i++ ){
         
         int r =(int)((float) (n-1-i)* rand()/(RAND_MAX+1.0)); 
         int temp = order_table[n-1-i]; //last 
         order_table[n-1-i] = order_table[r]; 
         order_table[r] = temp;
         }

    for (int i = 0; i < s; i++)
        order[i] = order_table[n-1-i];
    
    delete [] order_table;
    
    return;
    }

#endif
