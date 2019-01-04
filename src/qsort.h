// Init at Sept 9, 2016
// Bioinformatics Group, Single-Cell Research Center, QIBEBT, CAS
// Updated at June 9, 2018
// Updated by Gongchao Jing
// Quick sort for res
// Top K algorithm
#include <iostream>
#include <fstream>
using namespace std;

#ifndef QSORT_H
#define QSORT_H

//qsort for array ptr
int qsort_swap(float * score, int * cand, int a, int b){
    if (a == b) return 0;
    float tmp_f = score[a];
    score[a] = score[b];
    score[b] = tmp_f;
    
    if (cand == NULL) return 0;
    
    int tmp_i = cand[a];
    cand[a] = cand[b];
    cand[b] = tmp_i;
    return 0;
    }

int qsort_partition(float * score, int * cand, int begin, int end){
    
    float x = score[end -1];
    int i = begin - 1;
    for (int j = begin; j < end - 1; j ++)
        if (score[j] <= x){
            i ++;
            qsort_swap(score, cand, i, j);
            }
    qsort_swap(score, cand, i + 1, end - 1);
    
    return i + 1;
    }

int qsort_cand(float * score, int * cand, int begin, int end){
    if (begin < end - 1){
        int q = qsort_partition(score, cand, begin, end);
        qsort_cand(score, cand, begin, q);
        qsort_cand(score, cand, q + 1, end);
        }
    return 0;
    }
    
    
int insertK(float *score, int *cand, int cand_size, int swap_add)
{
	float tmp_f;int tmp_i; 
	for(int i=cand_size-1; i>=1;i--){
			if(score[swap_add]>=score[i-1]){
				tmp_f = score[i];tmp_i = cand[i];
				score[i]=score[swap_add]; cand[i] = cand[swap_add];
				score[swap_add] = tmp_f; cand[swap_add]=tmp_i;
				break;
				}
			else{
				tmp_f=score[i];tmp_i=cand[i];
				score[i]=score[i-1];cand[i]=cand[i-1];
				score[i-1]=score[swap_add];cand[i-1]=cand[swap_add];
				continue;
				}
		}
	
	}
	
int qsort(float * score, int * cand, int begin, int end,int cand_size){
    qsort_cand(score,cand,0,cand_size);
	for(int i=cand_size;i<end;i++){
		if(score[i]>score[cand_size-1])
			continue;
		else{
				insertK(score,cand,cand_size,i);
			}
		}
    return 0;
    }

#endif
