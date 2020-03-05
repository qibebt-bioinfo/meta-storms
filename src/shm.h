#include <iostream>

#include <sys/sem.h>
#include <sys/shm.h>
#include <sys/ipc.h>

#ifndef SHM_H
#define SHM_H

//debug

#define DBKEY_ST 4425
#define DBKEY_FH 4426
#define DBKEY_AG 4427

#define DBKEY_ST_FUNC 4435
#define DBKEY_FH_FUNC 4436
#define DBKEY_AG_FUNC 4437

#define DBKEY_ST_SP 4445
#define DBKEY_FH_SP 4446
#define DBKEY_AG_SP 4447

//MSE
/*
#define DBKEY_ST 1425
#define DBKEY_FH 1426
#define DBKEY_AG 1427

#define DBKEY_ST_FUNC 1435
#define DBKEY_FH_FUNC 1436
#define DBKEY_AG_FUNC 1437
*/
#define PATH_SIZE 1000

#define SHM_MODE (SHM_R | SHM_W)

typedef struct shm_args{
        
        char input_path[PATH_SIZE];
        char output_path[PATH_SIZE];
        short mode;
        short hit_n;
        short is_index;
        float min_s;
        }_shm_ag;

#endif
