// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at Dec 19, 2019
// Updated by Xiaoquan Su
// MetaDB_search

#include <sys/stat.h>
#include <omp.h>

#include "MetaDB.h"
#include "MetaDB_res.h"

#ifndef METADB_SEARCH_H
#define METADB_SEARCH_H

class _MetaDB_search : public _MetaDB{
      public:
             _MetaDB_search() : _MetaDB(){}
             _MetaDB_search(char db) : _MetaDB(db){
                                 MetaDB_Index = _MetaDB_Index(db);
                                 MetaDB_Comper = _MetaDB_Comper(db);
                                 }
             //Single sample
             vector <_MetaDB_Res> Index_Search(const char * queryfile, int n, int f, float t, float min_s, int mode, int coren);
             vector <_MetaDB_Res> Search(const char * queryfile, int n, float min_s, int mode, int coren);
    
             //Sample list
             vector <vector <_MetaDB_Res> > Index_Search(vector <string> queryfile, int n, int f, float t, float min_s, int mode, int coren);
             vector <vector <_MetaDB_Res> > Search(vector <string> queryfile, int n, float min_s, int mode, int coren);
    
             //Table format
             vector <vector <_MetaDB_Res> > Index_Search(_Table_Format * tablefile, int n, int f, float t, float min_s, int mode, int coren);
             vector <vector <_MetaDB_Res> > Search(_Table_Format * tablefile, int n, float min_s, int mode, int coren);  
             
             //Output the result in OTU table
             void Output_OTU_Table(const char * outfilename, const char * queryfile, string queryname, vector <_MetaDB_Res> res);                
	         void Output_Species_Table(const char * outfilename, const char * queryfile, string queryname, vector <_MetaDB_Res> res);	              
             
             
      private:
              int Load_Abd_Hdd(string id, float * abd);                          
              vector <_MetaDB_Res> Index_Search(float * q_abd, int n, int f, float t, float min_s, int mode, int coren);                      
              vector <_MetaDB_Res> Search(float * q_abd, int n, float min_s, int mode, int coren);                        

};

//Search releated methods
//Single sample
vector <_MetaDB_Res> _MetaDB_search::Index_Search(const char * queryfile, int n, int f, float t, float min_s, int mode, int coren){ //top n, index =  f * n
                    
                    float * q_abd = new float [MetaDB_Comper.Get_dim()];
                    MetaDB_Comper.Load_abd(queryfile, q_abd);
                    vector <_MetaDB_Res> res = Index_Search(q_abd, n, f, t, min_s, mode, coren);
                    delete [] q_abd;
                    return res; 
                    } 

vector <_MetaDB_Res> _MetaDB_search::Search(const char * queryfile, int n, float min_s, int mode, int coren){// top n
                    
                    float * q_abd = new float [MetaDB_Comper.Get_dim()];
                    MetaDB_Comper.Load_abd(queryfile, q_abd);
                    vector <_MetaDB_Res> res = Search(q_abd, n, min_s,mode, coren);
                    delete [] q_abd;
                    return res;
                    }

//Sample list
vector <vector <_MetaDB_Res> > _MetaDB_search::Index_Search(vector <string> queryfile, int n, int f, float t, float min_s, int mode, int coren){
    
    int sample_num = queryfile.size();

    vector <vector <_MetaDB_Res> > res_return;
    vector <_MetaDB_Res> * res = new vector <_MetaDB_Res> [sample_num];
    
    int single_sample_coren = 1;
    if (sample_num * 2 <= coren)
        single_sample_coren = coren / sample_num;
    
    omp_set_num_threads(coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < sample_num; i ++){
        res[i] = Index_Search(queryfile[i].c_str(), n, f, t, min_s, mode, single_sample_coren);
        }
    
    for (int i = 0; i < sample_num; i ++)
        res_return.push_back(res[i]);
    
    delete [] res;    
    return res_return;
}

vector <vector <_MetaDB_Res> > _MetaDB_search::Search(vector <string> queryfile, int n, float min_s, int mode, int coren){
    
    int sample_num = queryfile.size();
    
    vector <vector <_MetaDB_Res> > res_return;
    vector <_MetaDB_Res> * res = new vector <_MetaDB_Res> [sample_num];
    
    int single_sample_coren = 1;
    if (sample_num* 2 <= coren)
        single_sample_coren = coren / sample_num;
    
    omp_set_num_threads(coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < sample_num; i ++){
        
        res[i] = Search(queryfile[i].c_str(), n, min_s, mode, single_sample_coren);
        }
    
    for (int i = 0; i < sample_num; i ++)
        res_return.push_back(res[i]);
        
    delete [] res;
    return res_return;
    }

//Table format
vector <vector <_MetaDB_Res> > _MetaDB_search::Index_Search(_Table_Format * tablefile, int n, int f, float t, float min_s, int mode, int coren){
    
    int sample_num = tablefile->Get_Sample_Size();
       
    vector <vector <_MetaDB_Res> > res_return;
    vector <_MetaDB_Res> * res = new vector <_MetaDB_Res> [sample_num];
    
    int single_sample_coren = 1;
    if (sample_num * 2 <= coren)
        single_sample_coren = coren / sample_num;
    
    omp_set_num_threads(coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < sample_num; i ++){
           float * q_abd = new float [MetaDB_Comper.Get_dim()];
           MetaDB_Comper.Load_abd(tablefile, q_abd, i);
           res[i] = Index_Search(q_abd, n, f, t, min_s, mode, single_sample_coren);
           delete [] q_abd;
           }
    
    for (int i = 0; i < sample_num; i ++)
        res_return.push_back(res[i]);
    
    delete [] res;   
    return res_return;
    }

vector <vector <_MetaDB_Res> > _MetaDB_search::Search(_Table_Format * tablefile, int n, float min_s, int mode, int coren){
    
    int sample_num = tablefile->Get_Sample_Size();
    
    vector <vector <_MetaDB_Res> > res_return;
    vector <_MetaDB_Res> * res = new vector <_MetaDB_Res> [sample_num];
    
    int single_sample_coren = 1;
    if (sample_num * 2 <= coren)
        single_sample_coren = coren / sample_num;

    omp_set_num_threads(coren);
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < sample_num; i ++){
           float * q_abd = new float [MetaDB_Comper.Get_dim()];
           MetaDB_Comper.Load_abd(tablefile, q_abd, i);
           res[i] = Search(q_abd, n, min_s, mode, single_sample_coren);
           delete [] q_abd;
           }
    
    for (int i = 0; i < sample_num; i ++)
        res_return.push_back(res[i]);
    
    delete [] res;   
    return res_return;
    }

void _MetaDB_search::Output_OTU_Table(const char * outfilename, const char * queryfile, string queryname, vector <_MetaDB_Res> res){
                    
                    vector <string> otus;
                    for (int i = 0; i < MetaDB_Comper.Get_dim(); i ++)
                        otus.push_back("OTU_" + MetaDB_Comper.Get_Id(i));
                    
                    _Table_Format table(otus);
                    
                    float * q_abd = new float [MetaDB_Comper.Get_dim()];
                    MetaDB_Comper.Load_abd(queryfile, q_abd);
                    
                    vector <float> q_abd_vector(q_abd, q_abd + MetaDB_Comper.Get_dim());
                    
                    for (int i = 0; i < q_abd_vector.size(); i ++)
                        q_abd_vector[i] /= 100.0;
                        
                    table.Add_Abd(q_abd_vector, queryname);
                    
                    for (int i = 0; i < res.size(); i ++){
                        vector <float> res_abd_vector;
                        if (Mem_mode == 0){ //RAM                        
                           res_abd_vector = vector <float> ((res[i].Get_Sample())->Get_Abd(), (res[i].Get_Sample())->Get_Abd() + MetaDB_Comper.Get_dim());
                           }
                        else{ //HDD 
                              float * res_abd = new float [MetaDB_Comper.Get_dim()];
                              Load_Abd_Hdd((res[i].Get_Sample())->Get_ID(), res_abd);                   
                              res_abd_vector = vector <float> (res_abd, res_abd + MetaDB_Comper.Get_dim());
                             }
                              
                        for (int j = 0; j < res_abd_vector.size(); j ++)
                        res_abd_vector[j] /= 100.0;
                        table.Add_Abd(res_abd_vector, res[i].Get_Res_Sample_ID());
                    }
                    table.Filter_Empty();
                    table.Output_Table(outfilename);
                    
                    delete [] q_abd;
     
     }      
void _MetaDB_search::Output_Species_Table(const char * outfilename, const char * queryfile, string queryname, vector <_MetaDB_Res> res){

                    vector <string> otus;
                    for (int i = 0; i < MetaDB_Comper.Get_dim(); i ++)
                        otus.push_back(MetaDB_Comper.Get_Id(i));

                    _Table_Format table(otus);

                    float * q_abd = new float [MetaDB_Comper.Get_dim()];
                    MetaDB_Comper.Load_abd(queryfile, q_abd);

                    vector <float> q_abd_vector(q_abd, q_abd + MetaDB_Comper.Get_dim());

                    for (int i = 0; i < q_abd_vector.size(); i ++)
                        q_abd_vector[i] /= 100.0;

                    table.Add_Abd(q_abd_vector, queryname);

                    for (int i = 0; i < res.size(); i ++){
                        vector <float> res_abd_vector;
                        if (Mem_mode == 0){ //RAM                        
                           res_abd_vector = vector <float> ((res[i].Get_Sample())->Get_Abd(), (res[i].Get_Sample())->Get_Abd() + MetaDB_Comper.Get_dim());
                           }
                        else{ //HDD 
                              float * res_abd = new float [MetaDB_Comper.Get_dim()];
                              Load_Abd_Hdd((res[i].Get_Sample())->Get_ID(), res_abd);
                              res_abd_vector = vector <float> (res_abd, res_abd + MetaDB_Comper.Get_dim());
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
vector <_MetaDB_Res> _MetaDB_search::Index_Search(float * q_abd, int n, int f, float t, float min_s, int mode, int coren){ //top n, index =  f * n                                        
                    if (min_s < 0.01) min_s = 0.01; 
                    if (n <= 0) n = 1;                    
                    if (n > Samples.size()) n = Samples.size();
                    int cands_size = f * n;
                    if (cands_size > Samples.size()) cands_size = Samples.size();
                    
                    //q_abd count
                    int q_abd_count = 0;
                    for (int i = 0; i < MetaDB_Comper.Get_dim(); i ++)
                        q_abd_count += (q_abd[i] > 0) ? 1 : 0;
                    
                    
                    float * q_index_abd = new float [MetaDB_Index.Get_IndexN()];
                    MetaDB_Index.Make_Index(q_abd, q_index_abd);
                    
                    //Index calc
                    int * cands = new int [Samples.size()];
                    float * index_res = new float [Samples.size()];
    
                    omp_set_num_threads(coren);

                    #pragma omp parallel for schedule(dynamic, 1)
                    for (int i = 0; i < Samples.size(); i ++){
                        cands[i] = i;
                        index_res[i] = MetaDB_Index.Calc_Index(q_index_abd, Samples[i].Get_Index_Abd(), t);
                    }
    
                    //index sort
                    qsort(index_res, cands, 0, Samples.size(), cands_size);
    
                    //Index calc end
                                        
                    //Res calc
                    _MetaDB_Res * res = new _MetaDB_Res[cands_size];
                    
                    omp_set_num_threads(coren);
                   
                    #pragma omp parallel for schedule(dynamic, 1)
                    for (int i = 0; i < cands_size; i ++){
                        float sim = 0;
                        int res_abd_count = 0;
                        float * sample_abd = NULL;

                        if (Mem_mode == 0){ //RAM
                            sample_abd = Samples[cands[i]].Get_Abd();
                            res_abd_count = Samples[cands[i]].Get_OTU_Count();
                       		}
                        else { //HDD 
                             sample_abd =  new float [MetaDB_Comper.Get_dim()];
                             res_abd_count = Load_Abd_Hdd(Samples[cands[i]].Get_ID(), sample_abd);
                             }
                        sim = MetaDB_Comper.Calc_sim(q_abd, sample_abd, mode);
                        sim = (sim > 1.0) ? 1.0 : sim;
                        res[i] = _MetaDB_Res(&(Samples[cands[i]]), q_abd_count, res_abd_count, sim);     
                        
                        if (Mem_mode == 1) //HDD
                           delete [] sample_abd;                     
                        }
                    
                    _MetaDB_Res::Qsort(res, 0, cands_size);

                    vector <_MetaDB_Res> res_return;
                    for (int i = 0; i < n; i ++)
                        if (res[i].Sim_value >= min_s)
                           res_return.push_back(res[i]);
                    
                    delete [] q_index_abd;
                    delete [] cands;
                    delete [] index_res;
                    delete [] res;
                    
                    return res_return;
                    } 

vector <_MetaDB_Res> _MetaDB_search::Search(float * q_abd, int n, float min_s, int mode, int coren){// top n
                    if (min_s < 0.01) min_s = 0.01; 
                    if (n <= 0) n = 1;
                    if (n > Samples.size()) n = Samples.size();
                    
                    //q_abd count
                    int q_abd_count = 0;
                    for (int i = 0; i < MetaDB_Comper.Get_dim(); i ++)
                        q_abd_count += (q_abd[i] > 0) ? 1 : 0;
                    
                    _MetaDB_Res * res = new _MetaDB_Res[Samples.size()];
                    
                    omp_set_num_threads(coren);
                    
                    #pragma omp parallel for schedule(dynamic, 1)
                    for (int i = 0; i < Samples.size(); i ++){                                                
                        
                        float sim = 0;
                        int res_abd_count = 0;
                        float * sample_abd = NULL;

                        if (Mem_mode == 0){ //RAM
                            sample_abd = Samples[i].Get_Abd();
                            res_abd_count = Samples[i].Get_OTU_Count();
                        	}
                           
                        else { //HDD 
                             sample_abd =  new float [MetaDB_Comper.Get_dim()];
                             res_abd_count = Load_Abd_Hdd(Samples[i].Get_ID(), sample_abd);
                             }
                        sim = MetaDB_Comper.Calc_sim(q_abd, sample_abd, mode);
                        
                        sim = (sim > 1.0) ? 1.0 : sim;
                        res[i] = _MetaDB_Res(&(Samples[i]), q_abd_count, res_abd_count, sim);     
                        
                        if (Mem_mode == 1)//HDD
                           delete [] sample_abd;      
                         
                        }
                    
                    //drop sim = 0 matches, debug, for too many zeros 
                    float non_zero_cut = (min_s > 0.01) ? min_s : 0.01;
                    int non_zero_count = 0;
                    for (int i = 0; i < Samples.size(); i++)
                        if (res[i].Sim_value >= non_zero_cut) non_zero_count ++;
                    
                    _MetaDB_Res * non_zero_res = new _MetaDB_Res[non_zero_count];
                    
                    int non_zero_iterator = 0;
                    for (int i = 0; i < Samples.size(); i++)
                        if (res[i].Sim_value >= non_zero_cut){                                             
                                             non_zero_res[non_zero_iterator] = res[i];
                                             non_zero_iterator ++;
                                             }                    
                    if (non_zero_count < n) n = non_zero_count;
                    
                    //sort
                    /*
                    _MetaDB_Res::Qsort(res, 0, Samples.size());
    
                    vector <_MetaDB_Res> res_return;
                    for (int i = 0; i < n; i ++)
                    if (res[i].Sim_value >= min_s)
                           res_return.push_back(res[i]);                        
                    */
                    _MetaDB_Res::Qsort(non_zero_res, 0, non_zero_count);
    
                    vector <_MetaDB_Res> res_return;
                    for (int i = 0; i < n; i ++)
                           res_return.push_back(non_zero_res[i]); 
                    
                    delete [] res;
                    delete [] non_zero_res;
                    
                    return res_return;
                    }
                    
int _MetaDB_search:: Load_Abd_Hdd(string id, float * abd){
					if (Mem_mode == 0) return 0;
					memset(abd, 0, MetaDB_Comper.Get_dim() * sizeof(float));
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
