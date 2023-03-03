#ifndef _FPS_H_
#define _FPS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <sys/times.h> //these two h files are for linux
#include <unistd.h>

// mersenne twist
typedef struct fps_randgen fps_randgen;

struct fps_randgen
{
	unsigned mt[624];
	int 	 mti;
};


// Define a data structure for a literal in the SAT problem.
typedef struct lit lit;

struct lit {
    int clause_num;		//clause num, begin with 0
    int var_num;		//variable num, begin with 1
    int sense;			//is 1 for true literals, 0 for false literals.
};

typedef struct FPS FPS;

struct FPS{
    /*parameters of the instance*/
    int     num_vars;		//var index from 1 to num_vars 
    int     num_clauses;	//clause index from 0 to num_clauses-1
    int		max_clause_len;
    int		min_clause_len;
    int		formula_len;
    double	avg_clause_len;
    double 	ratio;
    int		ave_weight;   //significant score(sigscore) needed for aspiration
    int		delta_total_weight;
    int		threshold;
    float	p_scale;//w=w*p+ave_w*q
    float	q_scale;
    int		scale_ave;//scale_ave==ave_weight*q_scale
    int 	q_init;
    long long  mems_left;
    //cutoff
    int		max_tries;
    int		tries;
    int		max_flips;
    int		step;
    int* 	score_inc_vars;
    int* 	score_inc_flag;
    int 	score_inc_count;
    /* literal arrays */				
    lit**	var_lit;				//var_lit[i][j] means the j'th literal of var i.
    int*	var_lit_count;          //amount of literals of each var
    lit**	clause_lit;		        //clause_lit[i][j] means the j'th literal of clause i.
    int*	clause_lit_count; 	    // amount of literals in each clause		

    int*	cscc2;
    int   tabu_step;
    
    int   simple_cnum = 10;
    int   simple_vnum = 5;
    
    int   fps_cnum = 10;
    int   fps_vnum = 50;
    int   selected_nums;
    int*   sel_cs;
    int*   selected; // need init
    int*   selected_during_flip;
    int*   score2;
    int*   scores;
    int*   vars2;
    int   goodvar_stack2_num;
    int*   goodvar_stack2;
    int*   best_vars;
    
    int bandit_cnum = 20;
    int backward_step = 21; 
    double gama = 0.9, alpha = 1.0;
    double max_clause_score = 1000;
    int current_index;
    int if_exceed;
    int local_times; // need init
    int pre_unsat_nb;
    
    int* selected_clauses;
    int* selected_times; // need init
    int* sampled_clauses;
    double* clause_score; // need init
    
    
    
    /* Information about the variables. */
    int*    score;		
    int*	time_stamp;
    int*	fix;
    int*	cscc;

    /* Information about the clauses */			
    int*     clause_weight;	
    int*     sat_count;		
    int*	 sat_var;

    //unsat clauses stack
    int*	unsat_stack;		 //store the unsat clause number
    int		unsat_stack_fill_pointer;
    int*	index_in_unsat_stack;//which position is a clause in the unsat_stack

    //variables in unsat clauses
    int*	unsatvar_stack;		
    int		unsatvar_stack_fill_pointer;
    int*	index_in_unsatvar_stack;
    int*	unsat_app_count;		//a varible appears in how many unsat clauses

    //configuration changed decreasing variables (score>0 and confchange=1)
    int*	goodvar_stack;		
    int		goodvar_stack_fill_pointer;
    int*	already_in_goodvar_stack;

    //unit clauses preprocess
    lit*	unitclause_queue;		
    int		unitclause_queue_beg_pointer;
    int     unitclause_queue_end_pointer;
    int*    clause_delete;

    /* Information about solution */
    char*    cur_soln;	//the current solution, with 1's for True variables, and 0's for False variables
    fps_randgen randgen;
	char*	best_soln;
	int 	best_cost;
	// confl_trans
	int*   conflict_ct;
    int    in_conflict_sz;
	int*   in_conflict;

    // preprocess
    int     fix_var_ct;
    int     del_cls_ct;
};


void 	init_FPS(FPS*);
void 	reinit_FPS(FPS*);
void 	confl_trans(FPS*);


bool    local_search(FPS*);
int     build_instance(FPS*,char *filename);
void    print_solution(FPS*);
int     verify_sol(FPS*);
void 	alloc_memory(FPS*);
void    free_memory(FPS*);
void    settings(FPS*, char*);
void    fps_merseene_init(FPS*, int);
void 	update_after_build(FPS*);


#endif

