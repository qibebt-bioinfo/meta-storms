// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Mar 21, 2019
// Updated by Xiaoquan Su
// MetaDB_search_func

#include <sys/stat.h>
#include <omp.h>

#include "MetaDB_func.h"
#include "MetaDB_res.h"

#ifndef METADB_SEARCH_H
#define METADB_SEARCH_H

class _MetaDB_search_func : public _MetaDB_func{
      public:
             //Single sample
             vector <_MetaDB_Res> Search(const char * queryfile, int n, float min_s, int coren);
    
             //Sample list
             vector <vector <_MetaDB_Res> > Search(vector <string> queryfile, int n, float min_s, int coren);
    
             //Table format
             vector <vector <_MetaDB_Res> > Search(_Table_Format * tablefile, int n, float min_s, int coren);  
             
             //Output the result in OTU table
             void Output_Gene_Table(const char * outfilename, const char * queryfile, string queryname, vector <_MetaDB_Res> res);                
                          
             
      private:
                                          
              int Load_Abd_Hdd(string id, float * abd);                    
              vector <_MetaDB_Res> Search(float * q_abd, int n, float min_s, int coren);                        

};

//Search releated methods
//Single sample
vector <_MetaDB_Res> _MetaDB_search_func::Search(const char * queryfile, int n, float min_s, int coren){// top n
                    
                    float * q_abd = new float [Comp_tree_func.Get_GeneN()];
                    Comp_tree_func.Load_Gene_Count(queryfile, q_abd);
                    vector <_MetaDB_Res> res = Search(q_abd, n, min_s, coren);
                    delete [] q_abd;
                    return res;
                    }

//Sample list
vector <vector <_MetaDB_Res> > _MetaDB_search_func::Search(vector <string> queryfile, int n, float min_s, int coren){
    
    int sample_num = queryfile.size();
    
    vector <vector <_MetaDB_Res> > res_return;
    vector <_MetaDB_Res> * res = new vector <_MetaDB_Res> [sample_num];
    
    int single_sample_coren = 1;
    if (sample_num* 2 <= coren)
        single_sample_coren = coren / sample_num;
    
    omp_set_num_threads(coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < sample_num; i ++){
        res[i] = Search(queryfile[i].c_str(), n, min_s, single_sample_coren);
        }
    
    for (int i = 0; i < sample_num; i ++)
        res_return.push_back(res[i]);
        
    delete [] res;
    return res_return;
    }

//Table format
vector <vector <_MetaDB_Res> > _MetaDB_search_func::Search(_Table_Format * tablefile, int n, float min_s, int coren){
    
    int sample_num = tablefile->Get_Sample_Size();
    
    vector <vector <_MetaDB_Res> > res_return;
    vector <_MetaDB_Res> * res = new vector <_MetaDB_Res> [sample_num];
    
    int single_sample_coren = 1;
    if (sample_num * 2 <= coren)
        single_sample_coren = coren / sample_num;

    omp_set_num_threads(coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < sample_num; i ++){
           float * q_abd = new float [Comp_tree_func.Get_GeneN()];
           Comp_tree_func.Load_Gene_Count(tablefile, q_abd, i);
           res[i] = Search(q_abd, n, min_s, single_sample_coren);
           delete [] q_abd;
           }
    
    for (int i = 0; i < sample_num; i ++)
        res_return.push_back(res[i]);
    
    delete [] res;   
    return res_return;
    }

void _MetaDB_search_func::Output_Gene_Table(const char * outfilename, const char * queryfile, string queryname, vector <_MetaDB_Res> res){
                    
                    vector <string> otus;
                    for (int i = 0; i < Comp_tree_func.Get_GeneN(); i ++)
                        otus.push_back(Comp_tree_func.Get_Id(i));
                    
                    _Table_Format table(otus);
                    
                    float * q_abd = new float [Comp_tree_func.Get_GeneN()];
                    Comp_tree_func.Load_Gene_Count(queryfile, q_abd);
                    
                    vector <float> q_abd_vector(q_abd, q_abd + Comp_tree_func.Get_GeneN());
                    
                    for (int i = 0; i < Comp_tree_func.Get_GeneN(); i ++)
                        q_abd_vector[i] /= 100.0;
                        
                    table.Add_Abd(q_abd_vector, queryname);
                    
                    for (int i = 0; i < res.size(); i ++){
                        vector <float> res_abd_vector;
                        if (Mem_mode == 0){ //RAM                        
                           res_abd_vector = vector <float> ((res[i].Get_Sample())->Get_Abd(), (res[i].Get_Sample())->Get_Abd() + Comp_tree_func.Get_GeneN());
                           }
                        else{ //HDD 
                              float * res_abd = new float [Comp_tree_func.Get_GeneN()];
                              Load_Abd_Hdd((res[i].Get_Sample())->Get_ID(), res_abd);                   
                              res_abd_vector = vector <float> (res_abd, res_abd + Comp_tree_func.Get_GeneN());
                             }
                              
                        for (int j = 0; j < res_abd_vector.size(); j ++)
                        res_abd_vector[j] /= 100.0;
                        table.Add_Abd(res_abd_vector, res[i].Get_Res_Sample_ID());
                    }
                    table.Filter_Empty();
                    table.Output_Table(outfilename);
                    
                    delete [] q_abd;
     
     }      


//private

vector <_MetaDB_Res> _MetaDB_search_func::Search(float * q_abd, int n, float min_s, int coren){// top n
                    if (n <= 0) n = 1;
                    if (n > Samples.size()) n = Samples.size();
                    
                    //q_abd count
                    int q_abd_count = 0;
                    for (int i = 0; i < Comp_tree_func.Get_GeneN(); i ++)
                        q_abd_count += (q_abd[i] > 0) ? 1 : 0;
                    
                    _MetaDB_Res * res = new _MetaDB_Res[Samples.size()];
                    
                    omp_set_num_threads(coren);
                    
                    #pragma omp parallel for schedule(dynamic, 1)
                    for (int i = 0; i < Samples.size(); i ++){
                        float sim = 0;
                        float * sample_abd = NULL;
                        int res_abd_count = 0;
                                                
                         if (Mem_mode == 0){ //RAM
                            sample_abd = Samples[i].Get_Abd();
                            res_abd_count = Samples[i].Get_OTU_Count();
                        	}
                           
                        else { //HDD 
                             sample_abd =  new float [Comp_tree_func.Get_GeneN()];
                             res_abd_count = Load_Abd_Hdd(Samples[i].Get_ID(), sample_abd);
                             }
                        
                        //sim = 1 - Comp_tree_func.Calc_Dist_Cos(q_abd, sample_abd);
                        sim = Comp_tree_func.Calc_sim(q_abd, sample_abd, 0);
                        
                        sim = (sim > 1.0) ? 1.0 : sim;
                        res[i] = _MetaDB_Res(&(Samples[i]), q_abd_count, res_abd_count, sim);     
                        
                        //if (Mem_mode == 1)//HDD
                        //   delete [] sample_abd;                  
                        }
                    
                    //sort
                    _MetaDB_Res::Qsort(res, 0, Samples.size());
    
                    vector <_MetaDB_Res> res_return;
                    for (int i = 0; i < n; i ++)
                    if (res[i].Sim_value >= min_s)
                           res_return.push_back(res[i]);                        
                    
                    delete [] res;
                    
                    return res_return;
                    }

int _MetaDB_search_func::Load_Abd_Hdd(string id, float * abd){
					if (Mem_mode == 0) return 0;
					memset(abd, 0, Comp_tree_func.Get_GeneN() * sizeof(float));
	 				string hdd_path = Get_Sample_Hdd_Path(id);
					 //fgets
                    FILE * fptr = fopen(hdd_path.c_str(), "r");
                    if (fptr == NULL){
                       cerr << "Error: Cannot open input file: " << hdd_path << endl;
                       return -1;
                       }
                 
                  char * str_buffer = new char [BUFF];
                  vector <string> file_buffer;
                  
                   while(fgets(str_buffer, BUFF-1, fptr) != NULL){
                                        
                                        string a_line = str_buffer;
                                        file_buffer.push_back(a_line);
                                        
                                        }
                   delete [] str_buffer;
                   fclose(fptr);
                   int count = 0;
                   for(int i = 0; i < file_buffer.size(); i ++){
                                       int a_id = 0;
                                       float a_abd = 0;                                       
                                       sscanf(file_buffer[i].c_str(), "%d%f", &a_id, &a_abd);									   
                                       abd[a_id] = a_abd;
                                       count ++;
                                       }                              
					return count;
}


#endif
