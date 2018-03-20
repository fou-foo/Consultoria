
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h> 
#include <ctype.h>
#include <ilcplex\cplex.h>
/*foo*/
#define NUMROWS 90
#define NUMCOLS 800
#define NUMNZ    200 
#define INFBOUND 100
typedef struct _inputdata {
	int    plants;       /* Number of plants        */
	int    supp;         /* Number of suppliers     */
	int    scen;         /* Number of scenarios     */
	double **trancost;   /* Transportation costs    */
	int    *cap;         /* Supplier capacities     */
	double *fixedcost;   /* Supplier fixed cost     */
	int    **demand;     /* Demand scenarios        */
	double **exch;       /* Exchange rate scenarios */
	double *prob;        /* Probability values      */
	double  max_dem;       /* Maximium scenario dem   */
} inputdata;

typedef struct _soldata {
	int	*y;				 /* 0-1 variable   */
	double exp_value;    /* Expected Value */
    double av_trans;     /* Transportation cost (average) */
	double av_fixed;	 /* Fixed cost (average) */
	int sol_cap;		 /* Capacity of solution */
	double **dual_cost;  /* Matrix with dual costs */
	double *exp_dual_cost; /* Vector of expected " */
	int *ybest;
	double best_value;
	double **best_shadow;
        int best_iter;
	clock_t itime;
	clock_t btime;
} soldata;

typedef struct _datacplex {
	char     probname[16];  /* Problem name is max 16 characters */
	int      numcols;
	int      numrows;
	int      objsen;
	double   obj[NUMCOLS];
	double   rhs[NUMROWS];
	char     sense[NUMROWS];
	int      matbeg[NUMCOLS];
	int      matcnt[NUMCOLS];
	int      matind[NUMNZ];
	double   matval[NUMNZ];
	double   lb[NUMCOLS];
	double   ub[NUMCOLS];
	char     ctype[NUMCOLS];
	int		 index[NUMROWS];
	double	 values[NUMROWS];
	int      solstat;
	double   objval;
	double   z[NUMCOLS];
	double   pi[NUMROWS];
	double   slack[NUMROWS];
	double   dj[NUMCOLS];

} datacplex;

typedef struct _movedata {
	int type;
	int flag;
	int insertion;
	int deletion;
	double value;
	int *iter_in;
	int *iter_out;
	int **iter_swap;
} movedata;


void   input_data  (char *filename, inputdata *pdata);
void setproblemdata (datacplex *pplex, inputdata pdata,int l);
void initialize(inputdata pdata, soldata *psol);
int select_min(inputdata pdata,soldata *psol);
int evaluate(CPXENVptr env, CPXLPptr *lp,int iter,inputdata pdata,
				datacplex *pplex,soldata *psol,_int64 hcode, 
				int *ncoded,_int64 *coded_sol,double *coded_value);
void search(inputdata pdata,soldata *psol, CPXENVptr env, 
			CPXLPptr *lp, datacplex *pplex,movedata *move);
void form_list(inputdata pdata, soldata *psol, int *ins_list,
			   int *del_list, int **swap_list, int iter, 
			   movedata *move, _int64 *coded_sol, int *ncoded);
void select_move(CPXENVptr env,CPXLPptr *lp,inputdata pdata,soldata *psol,
				 datacplex *pplex, int *ins_list,int *del_list,
				 int **swap_list, movedata *move,int ins_length,
				 int del_length, int swap_length, int iter, _int64
				 *coded_sol, double *coded_value,int *ncoded);
void execute_move(CPXENVptr env,CPXLPptr *lp,datacplex *pplex,
				  inputdata pdata,soldata *psol,movedata *move, int iter);
_int64 code_sol(soldata *psol, inputdata pdata);
int is_there(_int64 hcode,_int64 *coded_sol,int ind, int *ncoded);
int find_place(_int64 hcode, _int64 *code_sol, int *ncoded);
void insert(_int64 hcode, _int64 *coded_sol, int *ncoded, int index);
void insert_value(double objetivo, double *coded_value, int *ncoded, int index);
_int64 n_choose_k(int n,int k);

int    *ivector    (int nl, int nh);
int    **imatrix   (int nrl, int nrh, int ncl, int nch);
double *dvector    (int nl, int nh);
double *dhvector    (_int64 nl, _int64 nh);
double **dmatrix   (int nrl, int nrh, int ncl, int nch);
_int64 *Hvector	   (_int64 nl, _int64 nh);

CPXLPptr *CPXLPvector(int nl, int nh);

void   free_ivector(int *v, int nl);
void   free_dvector(double *v, int nl);
void   free_imatrix(int **m, int nrl, int nrh, int ncl);
void   free_dmatrix(double **m, int nrl, int nrh, int ncl);
void free_Hvector(_int64 *v, _int64 nl);
void   nrerror     (char *error_text);

#define BIG = 999999


void write_file(char *filename,int n,int nsols,int **comb_matrix);



/*long int **limatrix(int nrl, int nrh, int ncl, int nch);
void free_limatrix(long int **m, int nrl, int nrh, int ncl);*/

int main(int argc,char **argv)

{

	int i,j;
	int inc,del;
	
	int		 status;

	inputdata pdata;     /* structure with the data of the problem */
	soldata psol;        /* structure with the solution of problem */
	datacplex pplex;     /* structure with cplex arrays            */
    movedata move;		 /* structure with move information	       */

	CPXENVptr env;
	CPXLPptr *lp; 
	double *rel_ben;

	//char in[] = "C:\\Users\\fou-f\\Desktop\\ROLOCTS\\EXAMPLE1.TXT";
	//char in[] = "C:\\Users\\fou-f\\Desktop\\ROLOCTS\\E87113.TXT";
	//char in[] = "C:\\Users\\fou-f\\Desktop\\ROLOCTS\\E543.TXT";
	//char in[] = "C:\\Users\\fou-f\\Desktop\\ROLOCTS\\0014010.TXT";
	/*	if (argc != 2) 
	{
		printf("usage: rolocts <input file>\n");
		exit(1);
	}*/

	input_data(argv[1], &pdata);

	
	/* Initialize solution structure */
	psol.y = ivector(1,pdata.supp);
	psol.exp_value=0;
	psol.sol_cap = 0;
	psol.dual_cost = dmatrix(1,pdata.scen,1,pdata.supp);
    psol.exp_dual_cost = dvector(1,pdata.supp); 
	psol.ybest = ivector(1,pdata.supp);
	psol.best_shadow = dmatrix(1,pdata.scen,1,pdata.supp);
	psol.best_value = 999999.0;
    rel_ben = dvector(1,pdata.supp);
	inc = 1;
	del = 0;

	/* Initialize move strucuture */
	move.type=0;
	move.flag = 0;
	move.insertion=0;
	move.deletion=0;
	move.value =0;
	move.iter_in = ivector(1,pdata.supp);
	move.iter_out = ivector(1,pdata.supp);
	move.iter_swap = imatrix(1,pdata.supp,1,pdata.supp);

    for(i=1; i<=pdata.supp; i++) 
	{
		move.iter_in[i] = 0;
		move.iter_out[i] = 0;
		for(j=1;j<=pdata.supp;j++) 
			move.iter_swap[i][j] = 0;
	}



	/* Initialize the CPLEX environment */
	env = CPXopenCPLEX (&status);
	if ( env == NULL )
	{
		char  errmsg[1024];
		fprintf (stderr, "Could not open CPLEX environment.\n");
		CPXgeterrorstring (env, status, errmsg);
		fprintf (stderr, "%s", errmsg);
	}
	lp = CPXLPvector(1, pdata.scen);
    for (i = 1; i <= pdata.scen; ++i)	
		lp[i] = NULL; 
   

   for (i=1; i<= pdata.scen; ++i) 
   {   
   	    setproblemdata (&pplex,pdata,i);
		lp[i] = CPXcreateprob(env, &status, pplex.probname);	
		status = CPXcopylp (env, lp[i], pplex.numcols, pplex.numrows, pplex.objsen, pplex.obj, pplex.rhs, pplex.sense, pplex.matbeg, pplex.matcnt, pplex.matind, pplex.matval,
                       pplex.lb, pplex.ub, NULL);
   }
	
   /* find a first feasible solution */	
   initialize(pdata,&psol); //solucion inicial de ceros y calculo de 'D'

   /* search for best solution */
   search(pdata,&psol,env,lp,&pplex,&move); 
return 0;
/* End Main */
}

/************************************************************
 *
 *     FIND PLACE OF AN ELEMENT IN THE SORTED ARRAY
 *
 ************************************************************/

 int find_place(_int64 hcode, _int64 *code_sol, int *ncoded)
 {
	 int i,nlow,nhigh;
	 _int64 key;
     
	 key = hcode;
	 nlow = 1;
	 nhigh = ncoded[0];

	 if(key < code_sol[nlow]) return nlow;
	 if(key > code_sol[nhigh]) return nhigh+1;

	 i = (nlow + nhigh)/2;
	 while(nlow < nhigh) {
		 if(key < code_sol[i]) {
			 nhigh = i;
		 }
		 else {
			 nlow = i;
		 }

		 if(nlow+1 == nhigh) return nhigh;

		 i = (nlow+nhigh)/2;
	 }
 }


 /***********************************************************
  *
  *         INSERT ELEMENT IN PROPER PLACE
  *
  ***********************************************************/

 void insert(_int64 hcode, _int64 *coded_sol, int *ncoded, int index)
 {
	 int i;
	 

	 if(index > ncoded[0]) {
		 coded_sol[index] = hcode;
	 }
	 else {
		 for(i=ncoded[0];i>=index;i--) {
			 coded_sol[i+1] = coded_sol[i];
		}

		coded_sol[index] = hcode;
	 }

 }

/**************************************************************
 *
 *           INSERT OBJECTIVE VALUE IN ARRAY
 *
 **************************************************************/

 void ins_value(double objetivo, double *coded_value, int *ncoded, int index)
 {
	 int i;
	 
	 

	 if(index > ncoded[0]) {
		 coded_value[index] = objetivo;
	 }
	 else {
		 for(i=ncoded[0];i>=index;i--) {
			 coded_value[i+1] = coded_value[i];
		 }
		 coded_value[index] = objetivo;
	 }
 }


/**********************************************************
 *
 *    SEARCH FOR SOLUTION IN DATA BASE
 *
 **********************************************************/
int is_there(_int64 hcode,_int64 *coded_sol,int ind, int *ncoded)

{
	_int64 key;
	int nlow, nhigh, i;

	ind = 0;
	key = hcode;
	nhigh = ncoded[0];
	nlow = 1;

	i = (nhigh + nlow)/2;
	while(nhigh - nlow > 1) {
		if(key == coded_sol[i]) {
			ind = i;
			return 1;
		}
		else {
			if(key < coded_sol[i]) {
				nhigh = i;
			}
			else {
				nlow = i;
			}
		}

		i = (nhigh + nlow)/2 ;
	}

	return 0;
}




/***********************************************************
 *
 *           CODE A GIVEN SOLUTION
 *
 **********************************************************/

 _int64 code_sol(soldata *psol, inputdata pdata)
 
 {
	 int i,n;
	 _int64 hcode;

	 hcode = 0;
	 n = pdata.supp;
	 for(i=1;i<=n;i++)  hcode += psol->y[i]*pow(2,n-i);

	 return hcode;
 }




/************************************************************
 *
 *             SEARCH PROCEDURE
 *
 ************************************************************/

 void search(inputdata pdata,soldata *psol, CPXENVptr env, CPXLPptr *lp, datacplex *pplex, movedata *move)
 {
	 int i;
	 int *ins_list,*del_list,**swap_list;
	 int iter,max_iter;
	 int ins_length,del_length,swap_length,increased;
     _int64 nsols;
	 _int64 *coded_sol;
	 double *coded_value;
	 _int64 hcode;
	 int *ncoded;
     double etime;
	 double ttime;
	 clock_t ftime;
	 
	 ///parametros IMPORTANTES
	 max_iter = 50;
	 iter = 1;
	
	 ins_length=3;
	 del_length=3;
	 swap_length= n_choose_k(pdata.supp,2)/4;
     nsols = pow(2,pdata.supp);
	 nsols = nsols/2;
	 
	 coded_sol = Hvector(1,nsols);
     coded_value = dvector(1,nsols);
	 ncoded = ivector(0,0);
	 ncoded[0] = 0;
	 hcode = 0;

	 ins_list = ivector(0,pdata.supp);
	 del_list = ivector(0,pdata.supp);
	 increased = 4*swap_length;
	 swap_list = imatrix(0,increased,1,2);

	 for(i=1; i<=pdata.supp; i++) 
	 {
		 if(psol->y[i]) move->iter_in[i] = 1;
	 }
	 evaluate(env, lp, iter, pdata, pplex, psol, hcode, ncoded, coded_sol, coded_value);
	 form_list(pdata,psol,ins_list,del_list,swap_list,iter,move, coded_sol,ncoded);
	 while(iter < max_iter) 
	 {
		 select_move(env,lp,pdata,psol,pplex,ins_list,del_list, swap_list,move,ins_length,del_length,swap_length, iter, coded_sol,coded_value,ncoded);
		 execute_move(env,lp,pplex,pdata,psol,move,iter);
		 form_list(pdata,psol,ins_list,del_list,swap_list,iter, move,coded_sol,ncoded);
			 iter +=1;
	 }
     ftime = clock();
	 etime = (double)(psol->btime - psol->itime)/CLOCKS_PER_SEC;
	 ttime = (double)(ftime - psol->itime)/CLOCKS_PER_SEC;
	 printf("  %4d  %7.2f  %7.2f  %10.6f\n", psol->best_iter,etime,ttime,psol->best_value);


	 /* end search procedure */
 }

 /*******************************************************
  *
  *                 EXECUTE MOVE
  *
  *******************************************************/

 void execute_move(CPXENVptr env,CPXLPptr *lp,datacplex *pplex,
		inputdata pdata,soldata *psol,movedata *move,int iter)

 {
	 int i, j;
	 int t_t_i, t_t_o, t_t_s;
	 int cnt, status, solstat1;
	 double objetivo;
	 t_t_i = pdata.supp/3;
	 t_t_o = pdata.supp/3;
	 t_t_s = n_choose_k(pdata.supp,2)/8; 
         
	 if(move->type ==1) 
	 {
		 psol->y[move->insertion] = 1;
		 psol->sol_cap += pdata.cap[move->insertion];
		 move->iter_out[move->insertion] = iter + t_t_i;
	 }
	 if(move->type == 2) 
	 {
		 psol->y[move->deletion] = 0;
		 psol->sol_cap -= pdata.cap[move->deletion];
		 move->iter_in[move->deletion] = iter + t_t_o;
	 }
	 if(move->type == 3) 
	 {
		 psol->y[move->insertion] = 1;
		 psol->y[move->deletion] = 0;
		 move->iter_swap[move->deletion][move->insertion] = iter + t_t_s;
		 move->iter_swap[move->insertion][move->deletion] = iter + t_t_s;
		 psol->sol_cap += pdata.cap[move->insertion]-pdata.cap[move->deletion];
	 }

	 if(move->flag == 0) 
	 {
		for(j=1;j<=pdata.scen;j++) 
		{
			 for(i=1;i<=pdata.supp;i++) 
			 {
				 psol->dual_cost[j][i] = psol->best_shadow[j][i];
			 }
		}
	 }
	 if(move->flag) 
	 {
		 for(i=1;i<=pdata.supp;i++) 
		 {
				pplex->index[i-1] = pdata.plants + i -1;
				pplex->values[i-1] = pdata.cap[i]*psol->y[i];
		 }
		 cnt = pdata.supp;
         for(j=1;j<=pdata.scen;j++)
		 {
				CPXchgrhs(env, lp[j], cnt, pplex->index, pplex->values);
		}
		 for(j=1;j<=pdata.scen; ++j) 
		 {
				status = CPXdualopt(env,lp[j]);
				if ( status ) 
				{
    				fprintf (stderr, "Failed to optimize LP dual.\n");
				}
			}
			for(j=1;j<=pdata.scen;j++)
			{
				status = CPXsolution (env, lp[j], &solstat1, &objetivo, pplex->z,  pplex->pi, pplex->slack, pplex->dj);
				if ( status ) 
				{
						fprintf (stderr, "Failed to obtain solution simplex.\n");
				}

				for(i=1;i<=pdata.supp;i++) 
				{
					psol->dual_cost[j][i] = pplex->pi[i];
				}
			}
	 }
	 if(move->value < psol->best_value)
	 {
		 psol->best_value = move->value;
		 psol->best_iter = iter;
		 psol->btime = clock();
	 }

	 if(move->type==1)
		 printf("Iteration %4d Insertion, supplier %3d, best value is %f\n", iter,move->insertion,psol->best_value);
	 if(move->type==2)
		 printf("Iteration %4d Deletion, supplier %3d, best value is %f\n", iter,move->deletion,psol->best_value);
	 if(move->type==3)
		 printf("Iteration %4d Swapp, deleted %3d, inserted %3d best value is %f\n", iter,move->deletion,move->insertion,psol->best_value);
	 for(i=1;i<=pdata.supp;i++) 
	 {
		 if(psol->y[i]) printf(" %2d",i);
	 }
	 printf(" \n");
 }




/********************************************************
 *
 *                SELECT MOVE
 *
 *******************************************************/

void select_move(CPXENVptr env,CPXLPptr *lp,inputdata pdata,soldata *psol,
				 datacplex *pplex, int *ins_list,int *del_list,
				 int **swap_list,movedata *move, int ins_length,
				 int del_length, int swap_length, int iter, _int64
				 *coded_sol, double *coded_value,int *ncoded)
{
	int i, j, k;
	int cnt, status, solstat1;
	double objetivo;
	double av_cost, bvar, pprob;
	double *trans, *fixed, *objvec;
	_int64 hcode;
	int ind, index;

	move->value = 1000000.0;
    trans = dvector(1,pdata.scen);
	fixed = dvector(1,pdata.scen);
	objvec = dvector(1,pdata.scen);
	ind = 0;

	/* Examine insertions */
	for(k=1;k<=ins_list[0];k++) 
	{ 	
		psol->y[ins_list[k]] = 1;
		hcode = code_sol(psol,pdata);
		if(!is_there(hcode,coded_sol,ind,ncoded))
		{
			for(i=1;i<=pdata.supp;i++) 
			{
				pplex->index[i-1] = pdata.plants + i -1;
				pplex->values[i-1] = pdata.cap[i]*psol->y[i];
			}
			cnt = pdata.supp;
			for(j=1;j<=pdata.scen;j++)
			{
				CPXchgrhs(env, lp[j], cnt, pplex->index, pplex->values);
			}	
			for(j=1; j<=pdata.scen; ++j) 
			{
				status = CPXdualopt(env,lp[j]);
				if ( status ) 
				{
    				fprintf (stderr, "Failed to optimize LP.\n");
				}
			}
			psol->exp_value = 0;
			for(j=1; j<=pdata.scen; j++) 
			{
				status = CPXsolution (env, lp[j], &solstat1, &objetivo, pplex->z,  pplex->pi, pplex->slack, pplex->dj);
				if ( status )
				{
						fprintf (stderr, "Failed to obtain solution simplex.\n");
				}

				for(i=1;i<=pdata.supp;i++) 
				{
					psol->dual_cost[j][i] = pplex->pi[i];
				}
				trans[j] = objetivo;
				fixed[j] = 0;
				for(i=1;i<=pdata.supp;i++)
				{
					fixed[j] += pdata.exch[i][j]*pdata.fixedcost[i]*psol->y[i];
				}
				objvec[j] = trans[j]+fixed[j];
			}
			psol->av_trans = 0;
			psol->av_fixed = 0;
			for(j=1;j<=pdata.scen;j++) 
			{
				psol->av_trans += pdata.prob[j]*trans[j];
				psol->av_fixed += pdata.prob[j]*fixed[j];
			}
			av_cost = psol->av_trans + psol->av_fixed;
    
			bvar = 0;
			pprob = 0;
			for(j=1;j<=pdata.scen;j++)
			{
				if(objvec[j] > av_cost) 
				{
					bvar += pdata.prob[j]*pow(objvec[j] - av_cost,2);
					pprob += pdata.prob[j];
				}
			}
			bvar = bvar / pprob;
			psol->exp_value = av_cost;
			objetivo = psol->exp_value;
			hcode = code_sol(psol,pdata);
			index = find_place(hcode,coded_sol,ncoded);
			insert(hcode,coded_sol,ncoded,index);
			ins_value(objetivo,coded_value,ncoded,index);
			ncoded[0] += 1;
			if(objetivo < move->value) 
			{
				if( (iter >= move->iter_in[ins_list[k]]) || (objetivo < psol->best_value) ) 
				{
					move->type=1;
					move->flag = 0;
					move->deletion = 0;
					move->insertion = ins_list[k];
					move->value = objetivo;
					for(j=1;j<=pdata.scen;j++) 
					{
						for(i=1;i<=pdata.supp;i++) 
						{
							psol->best_shadow[j][i] = psol->dual_cost[j][i];
						}
					}
				}
			}

			/* restore problem */
			psol->y[ins_list[k]] = 0;
		}
		else 
		{
			index = find_place(hcode,coded_sol,ncoded);
			objetivo = coded_value[index];
			if(objetivo < move->value) 
			{
				if( (iter >= move->iter_in[ins_list[k]]) || (objetivo < psol->best_value) )
				{
					move->type=1;
					move->flag = 1;
					move->deletion = 0;
					move->insertion = ins_list[k];
					move->value = objetivo;
				}
			}

			/* restore problem */
			psol->y[ins_list[k]] = 0;
		}
	}
	/* Examine deletions */
    for(k=1;k<=del_list[0];k++)
	{    
		psol->y[del_list[k]] = 0;
		hcode = code_sol(psol,pdata);
		if(!is_there(hcode,coded_sol,ind,ncoded))
		{
			for(i=1;i<=pdata.supp;i++) 
			{
				pplex->index[i-1] = pdata.plants + i -1;
				pplex->values[i-1] = pdata.cap[i]*psol->y[i];
			}
			cnt = pdata.supp;
			for(j=1;j<=pdata.scen;j++) 
			{
				CPXchgrhs(env, lp[j], cnt, pplex->index, pplex->values);
			}	
			for(j=1;j<=pdata.scen; ++j) 
			{
				status = CPXdualopt(env,lp[j]);
				if ( status )
				{
    				fprintf (stderr, "Failed to optimize LP.\n");
				}
			}
			psol->exp_value = 0;
			for(j=1;j<=pdata.scen;j++) 
			{
				status = CPXsolution (env, lp[j], &solstat1, &objetivo, pplex->z, pplex->pi, pplex->slack, pplex->dj);
				if ( status ) 
				{
						fprintf (stderr, "Failed to obtain solution.\n");
				}
				for(i=1;i<=pdata.supp;i++) 
				{
					psol->dual_cost[j][i] = pplex->pi[i];
				}
				trans[j] = objetivo;
				fixed[j] = 0;
				for(i=1;i<=pdata.supp;i++)
				{
					fixed[j] += pdata.exch[i][j]*pdata.fixedcost[i]*psol->y[i];
				}
				objvec[j] = trans[j]+fixed[j];
			}

			psol->av_trans = 0;
			psol->av_fixed = 0;
			for(j=1;j<=pdata.scen;j++)
			{
				psol->av_trans += pdata.prob[j]*trans[j];
				psol->av_fixed += pdata.prob[j]*fixed[j];
			}
			av_cost = psol->av_trans + psol->av_fixed;
			bvar = 0;
			pprob = 0;
			for(j=1;j<=pdata.scen;j++) 
			{
				if(objvec[j] > av_cost)
				{
					bvar += pdata.prob[j]*pow(objvec[j] - av_cost,2);
					pprob += pdata.prob[j];
				}
			}
			bvar = bvar / pprob;
			psol->exp_value = av_cost;
			objetivo = psol->exp_value;
			hcode = code_sol(psol,pdata);
			index = find_place(hcode,coded_sol,ncoded);
			insert(hcode,coded_sol,ncoded,index);
			ins_value(objetivo,coded_value,ncoded,index);
			ncoded[0] += 1;

			if(objetivo < move->value)
			{
				if( (iter >= move->iter_out[del_list[k]]) || (objetivo < psol->best_value) )
				{
					move->type=2;
					move->flag = 0;
					move->deletion = del_list[k];
					move->insertion = 0;
					move->value = objetivo;
					for(j=1;j<=pdata.scen;j++) 
					{
						for(i=1;i<=pdata.supp;i++) 
						{
							psol->best_shadow[j][i] = psol->dual_cost[j][i];
						}
					}
				}
			}

			/* restore problem */
			psol->y[del_list[k]]=1;
		}
		else 
		{
			index = find_place(hcode,coded_sol,ncoded);
			objetivo = coded_value[index];
			if(objetivo < move->value) 
			{
				if( (iter >= move->iter_out[del_list[k]]) || (objetivo < psol->best_value) ) {
					move->type=2;
					move->flag = 1;
					move->deletion = del_list[k];
					move->insertion = 0;
					move->value = objetivo;
				}
			}

			/* restore problem */
			psol->y[del_list[k]]=1;
		}
	}
	/*  Examine swapps */
	for(k=1;k<=swap_list[0][1];k++) 
	{    
		psol->y[swap_list[k][1]] = 0;
		psol->y[swap_list[k][2]] = 1;
		hcode = code_sol(psol,pdata);
		if(!is_there(hcode,coded_sol,ind,ncoded))
		{
			for(i=1;i<=pdata.supp;i++) 
			{
				pplex->index[i-1] = pdata.plants + i -1;
				pplex->values[i-1] = pdata.cap[i]*psol->y[i];
			}
			cnt = pdata.supp;
			for(j=1;j<=pdata.scen;j++) 
			{
				CPXchgrhs(env, lp[j], cnt, pplex->index, pplex->values);
			}
			for(j=1;j<=pdata.scen; ++j) 
			{
				status = CPXdualopt(env,lp[j]);
				if ( status ) 
				{
    				fprintf (stderr, "Failed to optimize LP dual .\n");
				}
			}

			psol->exp_value = 0;
			for(j=1;j<=pdata.scen;j++) 
			{
				status = CPXsolution (env, lp[j], &solstat1, &objetivo, pplex->z, pplex->pi, pplex->slack, pplex->dj);
				if ( status ) 
				{
						fprintf (stderr, "Failed to obtain solution simplex.\n");
				}
				for(i=1;i<=pdata.supp;i++) 
				{
					psol->dual_cost[j][i] = pplex->pi[i];
				}
				trans[j] = objetivo;
				fixed[j] = 0;
				for(i=1;i<=pdata.supp;i++) 
				{
					fixed[j] += pdata.exch[i][j]*pdata.fixedcost[i]*psol->y[i];
				}
				objvec[j] = trans[j]+fixed[j];
			}
			psol->av_trans = 0;
			psol->av_fixed = 0;
			for(j=1;j<=pdata.scen;j++) 
			{
				psol->av_trans += pdata.prob[j]*trans[j];
				psol->av_fixed += pdata.prob[j]*fixed[j];
			}
			av_cost = psol->av_trans + psol->av_fixed;
			bvar = 0;
			pprob = 0;
			for(j=1;j<=pdata.scen;j++) 
			{
				if(objvec[j] > av_cost) 
				{
					bvar += pdata.prob[j]*pow(objvec[j] - av_cost,2);
					pprob += pdata.prob[j];
				}
			}
			bvar = bvar / pprob;
			psol->exp_value = av_cost;
			objetivo = psol->exp_value;
			hcode = code_sol(psol,pdata);
			index = find_place(hcode,coded_sol,ncoded);
			insert(hcode,coded_sol,ncoded,index);
			ins_value(objetivo,coded_value,ncoded,index);
			ncoded[0] += 1;
			if(objetivo < move->value) 
			{
				if(iter >= move->iter_swap[swap_list[k][1]][swap_list[k][2]] || objetivo < psol->best_value  )
				{
					move->type=3;
					move->flag = 0;
					move->deletion = swap_list[k][1];
					move->insertion = swap_list[k][2];
					move->value = objetivo;
					for(j=1;j<=pdata.scen;j++) 
					{
						for(i=1;i<=pdata.supp;i++) 
						{
							psol->best_shadow[j][i] = psol->dual_cost[j][i];
						}
					}
				}
			}

			/* restore problem */
			psol->y[swap_list[k][1]] = 1;
			psol->y[swap_list[k][2]] = 0;
		}
		else 
		{
			index = find_place(hcode,coded_sol,ncoded);
			objetivo = coded_value[index];
			if(objetivo < move->value) 
			{
				if(iter >= move->iter_swap[swap_list[k][1]][swap_list[k][2]] || objetivo < psol->best_value  ) 
				{
					move->type=3;
					move->flag = 1;
					move->deletion = swap_list[k][1];
					move->insertion = swap_list[k][2];
					move->value = objetivo;
				}
			}
			/* restore problem */
			psol->y[swap_list[k][1]] = 1;
			psol->y[swap_list[k][2]] = 0;
		}
	}
	free_dvector(trans,1);
	free_dvector(fixed,1);
	free_dvector(objvec,1);
}

/*****************************************************
 *
 * COMPUTE COMBINATIONS OF K CHOSEN FROM SET N
 *
 ****************************************************/

_int64 n_choose_k(int n,int k)
{
	int i;
	_int64 numerator, denominator;
    _int64 nk;

	if(k<=n-k) {
		numerator = 1;
		for(i=n-k+1;i<=n;i++) {
			numerator *= i;
		}

		denominator = 1;
		for(i=2;i<=k;i++) {
			denominator *= i;
		}
	}

	else {
		numerator = 1;
		for(i=k+1;i<=n;i++) {
			numerator *= i;
		}

		denominator = 1;
		for(i=2;i<=n-k;i++) {
			denominator *= i;
		}
	}
	
	nk = numerator/denominator;

	return(nk);
}

 /****************************************************
  *
  *     FORM LIST
  *
  ***************************************************/

void form_list(inputdata pdata, soldata *psol, int *ins_list, int *del_list, int **swap_list, int iter, 
			   movedata *move, _int64 *coded_sol,int *ncoded)
{
	int i, j, k, ibest, jbest;
	double vbest;
	int *ins,*del,**swap;
	double *rel_ben;
	double **swap_ben;
	int ind;
	int max_ins_length;
	int max_del_length;
	int max_swap_length;
	_int64 hcode;
	max_ins_length = 3;
	max_del_length = 3;
	max_swap_length = n_choose_k(pdata.supp,2)/4;
	ind = 0;
	ins = ivector(1,pdata.supp);
    del = ivector(1,pdata.supp);
	swap = imatrix(1,pdata.supp,1,pdata.supp);
	rel_ben = dvector(1,pdata.supp);
	swap_ben = dmatrix(1,pdata.supp,1,pdata.supp);
	for(i=1;i<=pdata.supp;i++) 
	{
		psol->exp_dual_cost[i] = 0;
		for(j=1;j<=pdata.scen;j++) 
			psol->exp_dual_cost[i] += pdata.prob[j]*psol->dual_cost[j][i];
		rel_ben[i] = psol->exp_dual_cost[i] / pdata.fixedcost[i];
	}
	for(i=1; i<=pdata.supp; i++) 
	{
		for(j=1;j<=pdata.supp;j++) 
			swap_ben[i][j] = rel_ben[j] - rel_ben[i];
	}
	/* candidates for insertion */
	for(j=1;j<=pdata.supp;j++)	
		ins[j] = 0;
	for(i=1;i<= max_ins_length;i++)
	{
		jbest = 0;
		vbest = 99999;
		for(j=1;j<=pdata.supp;j++)
		{
			if(psol->y[j] == 0 && ins[j] == 0) 
			{
				if(iter>= move->iter_in[j])
				{
					if(rel_ben[j] < vbest) 
					{
						vbest = rel_ben[j];
						jbest = j;
					}
				}
			}
		}
		if(jbest) 
			ins[jbest] = 1;
		if(jbest) 
		{
			psol->y[jbest]=1;
			hcode = code_sol(psol,pdata);
			if(is_there(hcode,coded_sol,ind,ncoded)) 
				max_ins_length += 1;
			psol->y[jbest] = 0;
		}
	}
	/* candidates for deletion */
	for(j=1; j<=pdata.supp; j++)
		del[j] = 0;
	for(i=1; i<= max_del_length; i++)
	{
		jbest = 0;
		vbest = -99999;
		for(j=1; j<=pdata.supp; j++)
		{
			if(psol->sol_cap - pdata.cap[j] >= pdata.max_dem)
			{
				if(psol->y[j] == 1 && del[j] == 0) 
				{
					if(iter>=move->iter_out[j])
					{
						if(rel_ben[j] > vbest) 
						{
							vbest = rel_ben[j];
							jbest = j;
						}
					}
				}
			}
		}
		if(jbest) 
			del[jbest] = 1;
		if(jbest)
		{
			psol->y[jbest]=0;
			hcode = code_sol(psol,pdata);
			if(is_there(hcode,coded_sol,ind,ncoded)) 
				max_del_length += 1;
			psol->y[jbest] = 1;
		}
	}

	/* Candidates for swapping */
	for(i=1;i<=pdata.supp;i++)
	{
		for(j=1;j<=pdata.supp;j++) 
		{
			swap[i][j] = 0;
		}
	}
	for(k=1;k<= max_swap_length;k++)
	{
		ibest = 0;
		jbest = 0;
		vbest = 99999;
		for(i=1; i<=pdata.supp; i++)
		{
			if(psol->y[i]==1) 
			{
				for(j=1;j<=pdata.supp;j++) 
				{
					if(psol->sol_cap + pdata.cap[j] - pdata.cap[i] >= pdata.max_dem)
					{
						if(psol->y[j]==0 && swap[i][j]==0) 
						{
							if(iter>=move->iter_swap[i][j])
							{
								if(swap_ben[i][j] < vbest)
								{
									vbest = swap_ben[i][j];
									ibest = i;
									jbest = j;
								}
							}
						}
					}
				}
			}
		}
		if(ibest && jbest) 
			swap[ibest][jbest] = 1;
		if(ibest && jbest) 
		{
			psol->y[jbest]=1;
			psol->y[ibest]=0; 
			hcode = code_sol(psol,pdata);
			if(is_there(hcode,coded_sol,ind,ncoded)) 
				max_swap_length += 1;
			psol->y[jbest] = 0;
			psol->y[ibest] = 1;
			}
	}
	k = 0;
	for(j=1; j<=pdata.supp; j++)
	{
		if(ins[j]) 
		{
			k +=1;
			ins_list[k] = j;
		}
	}
	ins_list[0] = k;
	k = 0;
	for(j=1; j<=pdata.supp; j++) 
	{
		if(del[j]) 
		{
			k+=1;
			del_list[k] = j;
		}
	}
	del_list[0] = k;
	k = 0;
	for(i=1;i<=pdata.supp;i++) 
	{
		for(j=1;j<=pdata.supp;j++) 
		{
			if(swap[i][j])
			{
				k +=1;
				swap_list[k][1] = i;
				swap_list[k][2] = j;
			}
		}
	}
	swap_list[0][1] = k;
/* dispose of the arrays used in procedure */
	free_ivector(ins,1);
    free_ivector(del,1);
	free_imatrix(swap,1,pdata.supp,1);
	free_dvector(rel_ben,1);
	free_dmatrix(swap_ben,1,pdata.supp,1);
	/* end form_list */
}





/**********************************************************************
 *									                                  *
 *				  INPUT DATA     			                          *
 *									                                  *
 **********************************************************************/

void input_data(char *filename, inputdata *pdata)

{
	int  i, j;          /* counters                                   */
	char a;             /* used to scan the file and skip headers     */
	FILE *fp;           /* file pointer                               */
	errno_t err_open;
	int scen_dem;

	err_open = fopen_s(&fp, filename, "r");
	if (err_open != 0) nrerror("Unable to open input file.");

	/* read size of the problem after skipping headears */
	a = getc(fp);
	while(!isdigit(a)) a = getc(fp);
    ungetc(a, fp);
	fscanf_s(fp,"%d %d %d", &pdata->plants, &pdata->supp, &pdata->scen);

	/* allocate memory */
	pdata->cap       = ivector(1, pdata->supp);
	pdata->fixedcost = dvector(1, pdata->supp);
	pdata->prob      = dvector(1, pdata->scen);
	pdata->demand    = imatrix(1, pdata->plants, 1, pdata->scen);
	pdata->exch      = dmatrix(1, pdata->supp, 1, pdata->scen);
	pdata->trancost  = dmatrix(1, pdata->supp, 1, pdata->plants);

	/* read transportation costs (supplier, plant) */
	a = getc(fp);
    while(!isdigit(a)) a = getc(fp);
    ungetc(a, fp);
	for (i = 1; i <= pdata->supp; ++i)
		for (j = 1; j <= pdata->plants; ++j)
			fscanf_s(fp, "%lf", &pdata->trancost[i][j]);

	/* read supplier capacity and fixed costs */
	a = getc(fp);
    while(!isdigit(a)) a = getc(fp);
    ungetc(a, fp);
	for (i = 1; i <= pdata->supp; ++i)
		fscanf_s(fp, "%d %lf", &pdata->cap[i], &pdata->fixedcost[i]);

	/* read demand (plant, scenario) */
	a = getc(fp);
    while(!isdigit(a)) a = getc(fp);
    ungetc(a, fp);
	for (i = 1; i <= pdata->plants; ++i)
		for (j = 1; j <= pdata->scen; ++j)
			fscanf_s(fp, "%d", &pdata->demand[i][j]);

	/* read exchange rate (supplier, scenario) */
	a = getc(fp);
    while(!isdigit(a)) a = getc(fp);
    ungetc(a, fp);
	for (i = 1; i <= pdata->supp; ++i)
		for (j = 1; j <= pdata->scen; ++j)
			fscanf_s(fp, "%lf", &pdata->exch[i][j]);

	/* read probability values */
	a = getc(fp);
    while(!isdigit(a)) a = getc(fp);
    ungetc(a, fp);
	for (i = 1; i <= pdata->scen; ++i)
		fscanf_s(fp, "%lf", &pdata->prob[i]);
	fclose(fp);

	/* initialize maximum demand  */
	pdata->max_dem = 0;

	/* compute maximum demand */
	for(i=1;i<=pdata->scen;i++) {
		scen_dem = 0;
		for(j=1;j<=pdata->plants;j++) {
			scen_dem += pdata->demand[j][i];
		}
		if(scen_dem > pdata->max_dem) pdata->max_dem = scen_dem;
	}

}

/**************************************************************
 *
 *        SET PROBLEM DATA
 *
 **************************************************************/

void setproblemdata (datacplex *pplex, inputdata pdata,int l)
{   
   int i, j, k;  /* auxiliary indices  */
   	
   /* strcpy (probname, "lpsub"); */
   pplex->numcols   =  pdata.plants * pdata.supp ;
   pplex->numrows   = pdata.plants + pdata.supp;
   pplex->objsen = 1;   /* The problem is minimization */
 
   /* define matbeg */
   pplex->matbeg[0] = 0;

   for(i=1;i <= pplex->numcols -1;i++)
   {
		   pplex->matbeg[i] = pplex->matbeg[i-1] + 2;
   }
      
   /* define matcnt */
   for(i=0;i<= pplex->numcols -1;i++) 
   {
		   pplex->matcnt[i] = 2;
   }

   /* define matind */
   k = 0;
   for(i=1;i<=pdata.supp;++i) 
   {
	   for(j=1;j<=pdata.plants;++j) 
	   {
			   pplex->matind[k] = j-1;
			   pplex->matval[k] = 1.0;
			   k++;
			   pplex->matind[k] = pdata.plants + (i-1);
			   pplex->matval[k] = 1.0;
			   k++;
	   }
   }
   
   /* define bounds on the variables*/
   for(j=0;j<= pplex->numcols - 1; ++j)
   {
	   pplex->lb[j] = 0.0;
	   pplex->ub[j] = INFBOUND;
	   pplex->ctype[j] = 'C';
   }

   /* define the right hand side and inequalities direction*/
   for (i=0;i<= pdata.plants -1;i++) 
   {
	   pplex->sense[i] = 'E';
 	   pplex->rhs[i] = pdata.demand[i+1][l];
   }

   for (i=pdata.plants; i<= pplex->numrows -1; i++) 
   {
	   pplex->sense[i] = 'L';
	   pplex->rhs[i] = pdata.cap [i- pdata.plants +1];
   }

   /* define objective function coefficients  */
   k=0;
   for(i=0;i<=pdata.supp - 1;i++) 
   {
	   for(j=0;j<=pdata.plants - 1;j++) 
	   {
		   pplex->obj[k] = (pdata.exch[i+1][l])*(pdata.trancost[i+1][j+1]);
		   k++;
	   }
   }
}

/********************************************************
 *
 *    INITIALIZE SOLUTION
 *
 *******************************************************/

void initialize(inputdata pdata, soldata *psol)
{
	int i;
	int capacity;
    i = 0;
    capacity = 0;
	psol->itime = clock();
	for(i=1;i<=pdata.supp;i++) 
		psol->y[i]=0;
	
	while(capacity < pdata.max_dem) 
	{
		i = select_min(pdata,psol);
		capacity += (psol->y[i])*(pdata.cap[i]);
	}
	psol->sol_cap = capacity;
	for(i=1;i<=pdata.supp;i++) 
		psol->ybest[i]=psol->y[i];
}

/********************************
 * Select supplier with min cost
 ********************************/

int select_min(inputdata pdata,soldata *psol)
{
	int i,imin;
	double min;

	min = 100000;
	imin = 0;
	for(i=1;i<=pdata.supp;i++) {
		if(psol->y[i]==0) {
			if(pdata.fixedcost[i]/pdata.cap[i] < min) {
				min = pdata.fixedcost[i]/pdata.cap[i];
				imin = i;
			}
		}
	}

	psol->y[imin] = 1;
	return(imin);
}




/*********************
 * Evaluation function
 *********************/

int evaluate(CPXENVptr env, CPXLPptr *lp,int iter,inputdata pdata,
				datacplex *pplex,soldata *psol,_int64 hcode, 
				int *ncoded,_int64 *coded_sol,double *coded_value)
{
	int i, j, cnt;
	int status;
	double objetivo;
	int solstat1;
	double av_cost;
	double *trans,*fixed,*objvec;

	trans = dvector(1,pdata.scen);
	fixed = dvector(1,pdata.scen);
	objvec = dvector(1,pdata.scen);

	for(i=1; i<=pdata.supp; i++) 
	{
		pplex->index[i-1] = pdata.plants + i -1;
		pplex->values[i-1] = pdata.cap[i]*psol->y[i];
	}
	cnt = pdata.supp;    
	for(j=1;j<=pdata.scen;j++)
	{
			CPXchgrhs(env, lp[j], cnt, pplex->index, pplex->values);
	}
	if(iter==1) //primera iteracion 
	{
		for (j=1; j<= pdata.scen; ++j) 
		{	
   			status = CPXprimopt (env, lp[j]);
   			if ( status ) 
			{
    			fprintf (stderr, "Failed to optimize LP simplex.\n");
    			goto TERMINATE;
			}
		}
	}
	else 
	{
		for(j=1; j<=pdata.scen; ++j) 
		{
			status = CPXdualopt(env,lp[j]);
			if ( status ) 
			{
    			fprintf (stderr, "Failed to optimize LP dual.\n");
    			goto TERMINATE;
			}
		}
	}
	psol->exp_value = 0;
	for(j=1 ; j<=pdata.scen; j++) 
	{
		status = CPXsolution(env, lp[j], &solstat1, &objetivo, pplex->z, pplex->pi, pplex->slack, pplex->dj);
		if ( status ) 
		{
				fprintf (stderr, "Failed to obtain solution simplex.\n");
				goto TERMINATE;
		}
		for(i=1;i<=pdata.supp;i++) 
			psol->dual_cost[j][i] = pplex->pi[i];
		trans[j] = objetivo;
		fixed[j] = 0;
		for(i=1;i<=pdata.supp;i++) 
		{
			fixed[j] += pdata.exch[i][j]*pdata.fixedcost[i]*psol->y[i];
		}
		objvec[j] = trans[j]+fixed[j];
	}
	psol->av_trans = 0;
	psol->av_fixed = 0;
	for(j=1;j<=pdata.scen;j++) 
	{
		psol->av_trans += pdata.prob[j]*trans[j];
        psol->av_fixed += pdata.prob[j]*fixed[j];
	}
	av_cost = psol->av_trans + psol->av_fixed;
	psol->exp_value = av_cost;
	objetivo = psol->exp_value;

	hcode = code_sol(psol,pdata);
	coded_sol[1] = hcode;
	coded_value[1] = objetivo;
	ncoded[0] = 1;
	free_dvector(trans,1);
	free_dvector(fixed,1);
	free_dvector(objvec,1);
	return 0;
	TERMINATE:

   /* Free up the problem as allocated by CPXloadlp, if necessary */

   for (i=1; i<= pdata.scen; i++) 
   {
	   if ( lp[i] != NULL ) {
		  status = CPXfreeprob (env, &lp[i]);
		  if ( status ) {
			fprintf (stderr, "CPXfreeprob failed, error code %d.\n", status);
		  }
	   }
   }
}

/**************************************************************
 *
 *                 WRITE RESULTS TO A FILE
 *
 *************************************************************/

void write_file(char *filename,int n,int nsols,int **comb_matrix)

 {
   int j,k;       /* loop counter */

   FILE *f2;
   errno_t err_open;
   err_open = fopen_s(&f2, filename, "w");
   if (f2 == NULL) nrerror("Unable to open results file...");

   for (k=1;k<=nsols;k++) {
	   for(j=1;j<=n;j++) {
			fprintf(f2, "%4d\\%4d\n", k, comb_matrix[k][j]);
	   }
   }

   
}

/*****************************************************
 *
 *       MEMORY ALLOCATION PROCEDURES
 *
 ****************************************************/

int **imatrix(int nrl, int nrh, int ncl, int nch)

 {
     int i;
     int **m;

     m = (int **)malloc((unsigned) (nrh - nrl + 1)*sizeof(int*));
     if (!m) nrerror("allocation failure 1 in imatrix()");
     m -= nrl;
     for(i = nrl; i <= nrh ; i++)  {
        m[i] = (int*)malloc((unsigned) (nch - ncl + 1)*sizeof(int));
        if (!m[i]) nrerror("allocation failure 2 in imatrix()");
        m[i] -= ncl;
     }
     return m;
 }

 void free_imatrix(int **m, int nrl, int nrh, int ncl)

 {
     int i;

     for (i=nrh; i>=nrl; i--)
       free((char*) (m[i] + ncl));
       free ((char*) (m + nrl));
 }

 int *ivector(int nl, int nh)

 {
     int *v;
     v = (int *)malloc((unsigned)(nh-nl+1)*sizeof(int));
     if (!v) nrerror("Memory allocation failure in intvector");
       return v - nl;
 }
 
 
 void free_ivector(int *v, int nl)
 {
    free((char*) (v + nl));
 }

 double *dvector(int nl, int nh)

{
	double *v;

	v = (double *)malloc((unsigned)(nh-nl+1)*sizeof(double));
	if(!v)nrerror("allocation failure in fvector()");
	return v - nl;
}

double *dhvector(_int64 nl, _int64 nh)

{
	double *v;

	v = (double *)malloc((unsigned)(nh-nl+1)*sizeof(double));
	if(!v)nrerror("allocation failure in fvector()");
	return v - nl;
}

void free_dvector(double *v, int nl)

{
	free((char*) (v + nl));
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)

{
	int i;
	double **m;

	m = (double **)malloc((unsigned)(nrh - nrl +1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in fmatrix()");
	m -= nrl;
	for(i = nrl; i <= nrh ; i++) {
		m[i]= (double*)malloc((unsigned)(nch - ncl +1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in fmatrix()");
		m[i] -= ncl;
	}
	return m;
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl)

{
	int i;

	for (i = nrh; i >= nrl; i--)
		free((char*) (m[i] + ncl));
	free((char *)(m + nrl));
}


 
 void nrerror(char *error_text)
 {
   fprintf(stderr,"\n%s\n",error_text);
   fprintf(stderr,"...now exiting to system...");
   exit(1);
 }

 CPXLPptr *CPXLPvector(int nl, int nh)

{
	CPXLPptr *v;

	v = (CPXLPptr *)malloc((unsigned)(nh-nl+1)*sizeof(CPXLPptr));
	if(!v)nrerror("allocation failure in CPXLPvector()");
	return v - nl;
}


 _int64 *Hvector  (_int64 nl, _int64 nh)
 {
	_int64 *v;

	v = (_int64 *)malloc((unsigned)(nh-nl+1)*sizeof(_int64));
	if(!v)nrerror("allocation failure in fvector()");
	return v - nl;
 }


 void free_Hvector(_int64 *v, _int64 nl)

{
	free((char*) (v + nl));
}
