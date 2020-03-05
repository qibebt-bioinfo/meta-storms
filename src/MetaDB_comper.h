// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Dec 19, 2019
// Updated by Xiaoquan Su, Yufeng Zhang
// Automatic comper

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <map>
#include <set>
#include <omp.h>
#include <stdio.h>

#include "comp_dynamic.h"
#include "comp_sam_func.h"

#ifndef METADB_COMPER_H
#define METADB_COMPER_H

class _MetaDB_Comper{
      
      public:
              _MetaDB_Comper (){
                            Database = DEFAULT_DB;
                            Init();
                            }
             
             _MetaDB_Comper (char db){
                            Database = db;
                            Init();
                            }
                         
             int Load_abd(const char * infilename, float * abd);
             int Load_abd(_Table_Format * table, float * abd, int sample);
             
             int Get_dim();
             string Get_Id(int n);
             
             float Calc_sim(float * Abd_1, float * Abd_2);
             float Calc_sim(float * Abd_1, float * Abd_2, int mode);
      
      private:
              char Database; //G: otu greengenes; M: sp metaphlan; F: KO
              
              _Comp_Tree comp_tree_otu;
              _Comp_Tree_Dynamic comp_tree_sp;
              _Comp_Tree_Func comp_tree_func;
              
              void Init(){
                   comp_tree_otu = _Comp_Tree('G');
                   comp_tree_sp = _Comp_Tree_Dynamic('M');
                   comp_tree_func = _Comp_Tree_Func('G'); 
                   }      
      };

int _MetaDB_Comper::Load_abd(const char * infilename, float * abd){
    
    switch (Database){
           case 'G' : return comp_tree_otu.Load_abd(infilename, abd); break;
           case 'M' : return comp_tree_sp.Load_abd(infilename, abd); break;
           case 'F' : return comp_tree_func.Load_Gene_Count(infilename, abd); break;
           default: break;
           }
    return 0;
    }


int _MetaDB_Comper::Load_abd(_Table_Format * table, float * abd, int sample){
    
    switch (Database){
           case 'G' : return comp_tree_otu.Load_abd(table, abd, sample); break;
           case 'M' : return comp_tree_sp.Load_abd(table, abd, sample); break;
           case 'F' : return comp_tree_func.Load_Gene_Count(table, abd, sample); break;
           default: break;
           }
    return 0;
    }

int _MetaDB_Comper::Get_dim(){
    
    switch (Database){
           case 'G' : return comp_tree_otu.Get_LeafN(); break;
           case 'M' : return comp_tree_sp.Get_LeafN(); break;
           case 'F' : return comp_tree_func.Get_GeneN(); break;
           default: break;
           }
    return 0;
    
    }

string _MetaDB_Comper::Get_Id(int n){
       
       switch (Database){
           case 'G' : return comp_tree_otu.Get_Id(n); break;
           case 'M' : return comp_tree_sp.Get_Id(n); break;
           case 'F' : return comp_tree_func.Get_Id(n); break;
           default: break;
           }
    return "NULL";
       }

float _MetaDB_Comper::Calc_sim(float * Abd_1, float * Abd_2){
      
      switch (Database){
           case 'G' : return comp_tree_otu.Calc_sim(Abd_1, Abd_2); break;
           case 'M' : return comp_tree_sp.Calc_sim(Abd_1, Abd_2); break;
           case 'F' : return comp_tree_func.Calc_sim(Abd_1, Abd_2, 0); break;
           default: break;
           }
    return 0;      
      }

float _MetaDB_Comper::Calc_sim(float * Abd_1, float * Abd_2, int mode){
      
      switch (Database){
           case 'G' : return comp_tree_otu.Calc_sim(Abd_1, Abd_2, mode); break;
           case 'M' : return comp_tree_sp.Calc_sim(Abd_1, Abd_2, mode); break;
           case 'F' : return comp_tree_func.Calc_sim(Abd_1, Abd_2, mode); break;
           default: break;
           }
    return 0;      
      }

#endif
