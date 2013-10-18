#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "pll/src/pll.h"
#define GLOBAL_VARIABLES_DEFINITION
#include "nnisearch.h"
//#include "phylolib.h"
#include "fmemopen.h"

#if !defined WIN32 && !defined _WIN32 && !defined __WIN32__
#include <sys/resource.h>
#endif
/**
 verbose mode, determine how verbose should the screen be printed.
 */
enum VerboseMode {
	VB_QUIET, VB_MIN, VB_MED, VB_MAX, VB_DEBUG
};

/**
 verbose level on the screen
 */
extern int verbose_mode;

/* program options */
double TOL_LIKELIHOOD_PHYLOLIB;
int numSmoothTree;
int fast_eval;
int fivebran;
/* program options */

/*
int treeReadLenString(const char *buffer, pllInstance *tr, boolean readBranches,
		boolean readNodeLabels, boolean topologyOnly) {
	FILE *stream;

	stream = fmemopen((char*) buffer, strlen(buffer), "r");
	int lcount = treeReadLen(stream, tr, readBranches, readNodeLabels,
			topologyOnly);
	fclose(stream);
	return lcount;
}
*/

int cmp_nni(const void* nni1, const void* nni2) {
	NNIMOVE* myNNI1 = (NNIMOVE*) nni1;
	NNIMOVE* myNNI2 = (NNIMOVE*) nni2;
	if (myNNI1->likelihood > myNNI2->likelihood)
		return 1;
	else if (myNNI1->likelihood < myNNI2->likelihood)
		return -1;
	else
		return 0;
}

int compareDouble(const void * a, const void * b) {
	if (*(double*) a > *(double*) b)
		return 1;
	else if (*(double*) a < *(double*) b)
		return -1;
	else
		return 0;
}

NNIMOVE *getNNIList(pllInstance* tr) {
	static NNIMOVE* nniList;
	if (nniList == NULL) {
		nniList = (NNIMOVE*) malloc(2 * (tr->mxtips - 3) * sizeof(NNIMOVE));
		assert( nniList != NULL );
	}
	return nniList;
}

NNIMOVE *getNonConflictNNIList(pllInstance* tr) {
	static NNIMOVE* nonConfNNIList;
	if (nonConfNNIList == NULL) {
		nonConfNNIList = (NNIMOVE*) malloc((tr->mxtips - 3) * sizeof(NNIMOVE));
		assert( nonConfNNIList != NULL );
	}
	return nonConfNNIList;
}

double doNNISearch(pllInstance* tr, partitionList *pr, int* nni_count, double* deltaNNI) {
	/* save the initial tree likelihood, used for comparison */
	double initLH = tr->likelihood;

	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

	/* data structure to store the initial tree topology + branch length */
	topol* savedTree = setupTopol(tr->mxtips);
	saveTree(tr, savedTree, numBranches);

	/* Initialize the NNI list that holds information about 2n-6 NNI moves */
	NNIMOVE* nniList = getNNIList(tr);

	/* Now fill up the NNI list */
	nodeptr p = tr->start->back;
	nodeptr q = p->next;
	int numBran = 0; // number of visited internal branches during NNI evaluation
	int numPosNNI = 0; // number of positive NNI branhces found
	while (q != p) {
		evalNNIForSubtree(tr, pr, q->back, nniList, &numBran, &numPosNNI, initLH);
		q = q->next;
	}

	/* If no better NNI is found return -1.0 */
	if (numPosNNI == 0)
		return -1.0;

	int totalNNIs = numBran;
	/* Make sure all 2n-6 NNIs were evalauted */
	assert( totalNNIs == (2 * (tr->mxtips - 3)));
	/* Sort the NNI list ascendingly according to the log-likelihood */
	qsort(nniList, totalNNIs, sizeof(NNIMOVE), cmp_nni);

	/* Generate a list of independent positive NNI */
	NNIMOVE* inNNIs = getNonConflictNNIList(tr);

	/* The best NNI is the first to come to the list */
	inNNIs[0] = nniList[totalNNIs - 1];

	/* Subsequently add positive NNIs that are non-conflicting with the previous ones */
	int numInNNI = 1; // size of the existing non-conflicting NNIs in the list;
	for (int k = totalNNIs - 2; k > totalNNIs - numPosNNI - 1; k--) {
		int conflict = PLL_FALSE;
		/* Go through all the existing non-conflicting NNIs to check whether the next NNI will conflict with one of them */
		for (int j = 0; j < numInNNI; j++) {
			if (nniList[k].p->number == inNNIs[j].p->number
					|| nniList[k].p->number == inNNIs[j].p->back->number) {
				conflict = PLL_TRUE;
				break;
			}
			if (nniList[k].p->back->number == inNNIs[j].p->number
					|| nniList[k].p->back->number
							== inNNIs[j].p->back->number) {
				conflict = PLL_TRUE;
				break;
			}
		}
		if (conflict) {
			continue;
		} else {
			inNNIs[numInNNI] = nniList[k];
			numInNNI++;
		}
	}

	/* Applying all independent NNI moves */
	int numNNI = numInNNI;
	int MAXSTEPS = 100;
	int rollBack = PLL_FALSE;
	for (int step = 1; step <= MAXSTEPS; step++) {
		for (int i = 0; i < numNNI; i++) {
			/* First, do the topological change */
			doOneNNI(tr, pr, inNNIs[i].p, inNNIs[i].nniType, TOPO_ONLY);
			/* Then apply the new branch lengths */
			for (int j = 0; j < numBranches; j++) {
				inNNIs[i].p->z[j] = inNNIs[i].z0[j];
				inNNIs[i].p->back->z[j] = inNNIs[i].z0[j];
				inNNIs[i].p->next->z[j] = inNNIs[i].z1[j];
				inNNIs[i].p->next->back->z[j] = inNNIs[i].z1[j];
				inNNIs[i].p->next->next->z[j] = inNNIs[i].z2[j];
				inNNIs[i].p->next->next->back->z[j] = inNNIs[i].z2[j];
				inNNIs[i].p->back->next->z[j] = inNNIs[i].z3[j];
				inNNIs[i].p->back->next->back->z[j] = inNNIs[i].z3[j];
				inNNIs[i].p->back->next->next->z[j] = inNNIs[i].z4[j];
				inNNIs[i].p->back->next->next->back->z[j] = inNNIs[i].z4[j];
			}
			/* Update the partial likelihood */
			// TODO is it needed here? Calling newviewGeneric is very time consuming
		    if(numBranches > 1 && !tr->useRecom) {
		      pllNewviewGeneric(tr, pr, inNNIs[i].p, PLL_TRUE);
		      pllNewviewGeneric(tr, pr, inNNIs[i].p->back, PLL_TRUE);
		    }
		    else {
		      pllNewviewGeneric(tr, pr, inNNIs[i].p, PLL_FALSE);
		      pllNewviewGeneric(tr, pr, inNNIs[i].p->back, PLL_FALSE);
		    }
		}
		pllTreeEvaluate(tr, pr, 1);
		/* new tree likelihood should not be smaller the likelihood of the computed best NNI */
		if (tr->likelihood < inNNIs[0].likelihood) {
			if (numNNI == 1) {
				printf(
						"ERROR: new logl=%10.4f after applying only the best NNI < best NNI logl=%10.4f\n",
						tr->likelihood, inNNIs[0].likelihood);
				exit(1);
			}
			printf(
					"%d NNIs logl=%10.4f/best NNI logl= %10.4f. Roll back tree...\n",
					numNNI, tr->likelihood, inNNIs[0].likelihood);

			/* restore all branch lengths */
			if ( ! restoreTree(savedTree, tr, pr) ) {
				printf("ERROR: failed to roll back tree \n");
				exit(1);
			} else {
				rollBack = PLL_TRUE;
			}
			numNNI = 1;
		} else {
			if (rollBack) {
				printf("New logl:%10.4f\n", tr->likelihood);
			}
			break;
		}

	}
	*nni_count = numNNI;
	*deltaNNI = (tr->likelihood - initLH) / numNNI;
	return tr->likelihood;
}

/** @brief Optimize the length of a specific branch with variable number of Newton Raphson iterations

    Optimize the length of the branch connecting \a p and \a p->back
    for each partition (\a tr->numBranches) in library instance \a tr.

    @param tr
      The library instance

    @param pr
      Partition list

    @param p
      Endpoints of branch to be optimized

    @param newzpercycle Maximal number of Newton Raphson iterations

*/
void _update(pllInstance *tr, partitionList *pr, nodeptr p)
{
  nodeptr  q;
  int i;
  double   z[PLL_NUM_BRANCHES], z0[PLL_NUM_BRANCHES];
  int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

  q = p->back;

  for(i = 0; i < numBranches; i++)
    z0[i] = q->z[i];

  if(numBranches > 1)
    makenewzGeneric(tr, pr, p, q, z0, IQTREE_NEWZPERCYCLE, z, PLL_TRUE);
  else
    makenewzGeneric(tr, pr, p, q, z0, IQTREE_NEWZPERCYCLE, z, PLL_FALSE);

  for(i = 0; i < numBranches; i++)
  {
    if(!tr->partitionConverged[i])
    {
      if(PLL_ABS(z[i] - z0[i]) > PLL_DELTAZ)
      {
        tr->partitionSmoothed[i] = PLL_FALSE;
      }

      p->z[i] = q->z[i] = z[i];
    }
  }
}

double doOneNNI(pllInstance *tr, partitionList *pr, nodeptr p, int swap, int evalType) {
	nodeptr q;
	nodeptr tmp;

	q = p->back;
	assert(!isTip(q->number, tr->mxtips));
	assert(!isTip(p->number, tr->mxtips));
	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;

	if (swap == 1) {
		tmp = p->next->back;
		hookup(p->next, q->next->back, q->next->z, numBranches);
		hookup(q->next, tmp, tmp->z, numBranches);
	} else {
		tmp = p->next->next->back;
		hookup(p->next->next, q->next->back, q->next->z, numBranches);
		hookup(q->next, tmp, tmp->z, numBranches);
	}

	if (evalType == TOPO_ONLY) {
		return 0.0;
	} else if (evalType == ONE_BRAN_OPT) {
	    if(numBranches > 1 && !tr->useRecom) {
	      pllNewviewGeneric(tr, pr,p, PLL_TRUE);
	      pllNewviewGeneric(tr, pr,q, PLL_TRUE);

	    }
	    else {
	      pllNewviewGeneric(tr, pr,p, PLL_FALSE);
	      pllNewviewGeneric(tr, pr,q, PLL_FALSE);
	    }
		_update(tr, pr, p);
		pllEvaluateGeneric(tr, pr, p, PLL_FALSE, PLL_FALSE);
	} else { // 5 branches optimization
	    if(numBranches > 1 && !tr->useRecom) {
	      pllNewviewGeneric(tr, pr,p, PLL_TRUE);
	      pllNewviewGeneric(tr, pr,q, PLL_TRUE);
	    }
	    else {
	      pllNewviewGeneric(tr, pr,p, PLL_FALSE);
	      pllNewviewGeneric(tr, pr,q, PLL_FALSE);
	    }
	    _update(tr, pr, p);
		nodeptr r; // temporary node poiter
		// optimize 2 branches at node p
		r = p->next;
	    if(numBranches > 1 && !tr->useRecom)
	      pllNewviewGeneric(tr, pr,r, PLL_TRUE);
	    else
	      pllNewviewGeneric(tr, pr,r, PLL_FALSE);
	    _update(tr, pr, r);
		r = p->next->next;
	    if(numBranches > 1 && !tr->useRecom)
	      pllNewviewGeneric(tr, pr, r, PLL_TRUE);
	    else
	      pllNewviewGeneric(tr, pr, r, PLL_FALSE);
	    _update(tr, pr, r);
		// optimize 2 branches at node q
		r = q->next;
	    if(numBranches > 1 && !tr->useRecom)
	      pllNewviewGeneric(tr, pr, r, PLL_TRUE);
	    else
	      pllNewviewGeneric(tr, pr, r, PLL_FALSE);
	    _update(tr, pr, r);
		r = q->next->next;
	    if(numBranches > 1 && !tr->useRecom)
	      pllNewviewGeneric(tr, pr, r, PLL_TRUE);
	    else
	      pllNewviewGeneric(tr, pr, r, PLL_FALSE);
	    _update(tr, pr, r);
	    pllEvaluateGeneric(tr, pr, r, PLL_FALSE, PLL_FALSE);
	}
	return tr->likelihood;

}

int evalNNIForBran(pllInstance* tr, partitionList *pr, nodeptr p, NNIMOVE* nniList, int* numBran,
		double curLH) {
	nodeptr q = p->back;
	assert(!isTip(p->number, tr->mxtips));
	assert(!isTip(q->number, tr->mxtips));
	int numBranches = pr->perGeneBranchLengths ? pr->numberOfPartitions : 1;
	int betterNNI = 0;
	int i;
	NNIMOVE nni0; // dummy NNI to store backup information
	nni0.p = p;
	nni0.nniType = 0;
	nni0.likelihood = curLH;
	for (i = 0; i < numBranches; i++) {
		nni0.z0[i] = p->z[i];
		nni0.z1[i] = p->next->z[i];
		nni0.z2[i] = p->next->next->z[i];
		nni0.z3[i] = q->next->z[i];
		nni0.z4[i] = q->next->next->z[i];
	}

//    double multiLH = 0.0;
//    if (nnicut->doNNICut) {
//        // compute likelihood of the multifurcating tree
//        for (i = 0; i < tr->numBranches; i++) {
//            p->z[i] = 1.0;
//            p->back->z[i] = 1.0;
//        }
//        // now compute the log-likelihood
//        // This is actually time consuming!!!
//        newviewGeneric(tr, p, FALSE);
//        newviewGeneric(tr, p->back, FALSE);
//        evaluateGeneric(tr, p, FALSE);
//        multiLH = tr->likelihood;
//        //printf("curLH - multiLH : %f - %f \n", curLH, multiLH);
//        for (i = 0; i < tr->numBranches; i++) {
//            p->z[i] = z0[i];
//            p->back->z[i] = z0[i];
//        }
//        // If the log-likelihood of the zero branch configuration is some delta log-likelihood smaller than
//        // the log-likelihood of the current tree then it is very likely that there no better NNI tree
//        if (curLH - multiLH > nnicut->delta_min)
//            return nni0;
//    }

	/* TODO Save the likelihood vector at node p and q */
	//saveLHVector(p, q, p_lhsave, q_lhsave);
	/* Save the scaling factor */

	/* do an NNI move of type 1 */
	double lh1;
	if (!fast_eval && !fivebran) {
		lh1 = doOneNNI(tr, pr, p, 1, ONE_BRAN_OPT);
	} else if (fivebran) {
		lh1 = doOneNNI(tr, pr, p, 1, FIVE_BRAN_OPT);
	} else {
		lh1 = doOneNNI(tr, pr, p, 1, NO_BRAN_OPT);
	}
	NNIMOVE nni1;
	nni1.p = p;
	nni1.nniType = 1;
	// Store the optimized branch lengths
	for (i = 0; i < numBranches; i++) {
		nni1.z0[i] = p->z[i];
		nni1.z1[i] = p->next->z[i];
		nni1.z2[i] = p->next->next->z[i];
		nni1.z3[i] = q->next->z[i];
		nni1.z4[i] = q->next->next->z[i];
	}
	nni1.likelihood = lh1;

	nniList[*numBran] = nni1;

	if (nni1.likelihood > curLH + TOL_LIKELIHOOD_PHYLOLIB) {
		betterNNI = 1;
	}

	/* Restore previous NNI move */
	doOneNNI(tr, pr, p, 1, TOPO_ONLY);
	/* Restore the old branch length */
	for (i = 0; i < numBranches; i++) {
		p->z[i] = nni0.z0[i];
		q->z[i] = nni0.z0[i];
		p->next->z[i] = nni0.z1[i];
		p->next->back->z[i] = nni0.z1[i];
		p->next->next->z[i] = nni0.z2[i];
		p->next->next->back->z[i] = nni0.z2[i];
		q->next->z[i] = nni0.z3[i];
		q->next->back->z[i] = nni0.z3[i];
		q->next->next->z[i] = nni0.z4[i];
		q->next->next->back->z[i] = nni0.z4[i];
	}

	/* do an NNI move of type 2 */
	double lh2;
	if (!fast_eval && !fivebran) {
		lh2 = doOneNNI(tr, pr, p, 2, ONE_BRAN_OPT);
	} else if (fivebran) {
		lh2 = doOneNNI(tr, pr, p, 2, FIVE_BRAN_OPT);
	} else {
		lh2 = doOneNNI(tr, pr, p, 2, NO_BRAN_OPT);
	}
	// Create the nniMove struct to store this move
	NNIMOVE nni2;
	nni2.p = p;
	nni2.nniType = 2;
	// Store the optimized and unoptimized central branch length
	for (i = 0; i < numBranches; i++) {
		nni2.z0[i] = p->z[i];
		nni2.z1[i] = p->next->z[i];
		nni2.z2[i] = p->next->next->z[i];
		nni2.z3[i] = q->next->z[i];
		nni2.z4[i] = q->next->next->z[i];
	}
	nni2.likelihood = lh2;

	nniList[*numBran + 1] = nni2;

	if (nni2.likelihood > curLH + TOL_LIKELIHOOD_PHYLOLIB) {
		betterNNI = 1;
	}

	/* Restore previous NNI move */
	doOneNNI(tr, pr, p, 2, TOPO_ONLY);
	/* Restore the old branch length */
	for (i = 0; i < numBranches; i++) {
		p->z[i] = nni0.z0[i];
		q->z[i] = nni0.z0[i];
		p->next->z[i] = nni0.z1[i];
		p->next->back->z[i] = nni0.z1[i];
		p->next->next->z[i] = nni0.z2[i];
		p->next->next->back->z[i] = nni0.z2[i];
		q->next->z[i] = nni0.z3[i];
		q->next->back->z[i] = nni0.z3[i];
		q->next->next->z[i] = nni0.z4[i];
		q->next->next->back->z[i] = nni0.z4[i];
	}
	/* TODO: this should be replace replace by a function which restores the partial likelihood */
    if(numBranches > 1 && !tr->useRecom) {
      pllNewviewGeneric(tr, pr,p, PLL_TRUE);
      pllNewviewGeneric(tr, pr,p->back, PLL_TRUE);
    }
    else {
      pllNewviewGeneric(tr, pr,p, PLL_FALSE);
      pllNewviewGeneric(tr, pr,p->back, PLL_TRUE);

    }
//    if (nnicut->doNNICut && (nnicut->num_delta) < MAX_NUM_DELTA) {
//        double lh_array[3];
//        lh_array[0] = curLH;
//        lh_array[1] = lh1;
//        lh_array[2] = lh2;
//
//        qsort(lh_array, 3, sizeof (double), compareDouble);
//
//        double deltaLH = lh_array[1] - multiLH;
//        if (deltaLH > nnicut->delta_min) {
//            printf("%f %f %f %f\n", multiLH, curLH, lh1, lh2);
//        }
//        nnicut->delta[nnicut->num_delta] = deltaLH;
//        (nnicut->num_delta)++;
//    }

	return betterNNI;
	/* Restore the likelihood vector */
	//restoreLHVector(p,q, p_lhsave, q_lhsave);
}

/*int test_nni(int argc, char * argv[])
 {
 struct rusage usage;
 struct timeval start, end;
 struct timeval start_tmp, end_tmp;
 getrusage(RUSAGE_SELF, &usage);
 start = usage.ru_utime;
 tree        * tr;
 if (argc < 2)
 {
 fprintf (stderr, "syntax: %s [binary-alignment-file] [prefix] \n", argv[0]);
 return (1);
 }
 tr = (tree *)malloc(sizeof(tree));

 read the binary input, setup tree, initialize model with alignment
 read_msa(tr,argv[1]);

 tr->randomNumberSeed = 665;

 Create random tree
 //makeRandomTree(tr);
 //printf("RANDOM TREE: Number of taxa: %d\n", tr->mxtips);
 //printf("RANDOM TREE: Number of partitions: %d\n", tr->NumberOfModels);

 printf("Creating parsimony tree ...\n");
 getrusage(RUSAGE_SELF, &usage);
 start_tmp = usage.ru_utime;
 makeParsimonyTree(tr);
 getrusage(RUSAGE_SELF, &usage);
 end_tmp = usage.ru_utime;
 printf("Parsimony tree created in %f seconds \n", (double) end_tmp.tv_sec + end_tmp.tv_usec/1.0e6 - start_tmp.tv_sec - start_tmp.tv_usec/1.0e6);

 //printf("PARSIMONY TREE: Number of taxa: %d\n", tr->mxtips);
 //printf("PARSIMONY TREE: Number of partitions: %d\n", tr->NumberOfModels);

 compute the LH of the full tree
 //printf ("Virtual root: %d\n", tr->start->number);

 int printBranchLengths=TRUE;
 Tree2String(tr->tree_string, tr, tr->start->back, printBranchLengths, TRUE, 0, 0, 0, SUMMARIZE_LH, 0,0);
 //fprintf(stderr, "%s\n", tr->tree_string);
 char parTree[100];
 strcpy (parTree, argv[1]);
 strcat(parTree,".parsimony_tree");
 FILE *parTreeFile = fopen(parTree, "w");
 fprintf(parTreeFile, "%s\n", tr->tree_string);
 fclose(parTreeFile);

 Model optimization
 printf("Optimizing model parameters and branch lengths ... \n");
 getrusage(RUSAGE_SELF, &usage);
 start_tmp = usage.ru_utime;
 modOpt(tr, 0.1);
 getrusage(RUSAGE_SELF, &usage);
 end_tmp = usage.ru_utime;
 printf("Model parameters and branch lengths optimized in %f seconds \n", (double) end_tmp.tv_sec + end_tmp.tv_usec/1.0e6 - start_tmp.tv_sec - start_tmp.tv_usec/1.0e6);
 evaluateGeneric(tr, tr->start, FALSE);
 printf("Likelihood of the parsimony tree: %f\n", tr->likelihood);

 8 rounds of branch length optimization
 //smoothTree(tr, 32);
 //evaluateGeneric(tr, tr->start, TRUE);
 //printf("Likelihood after branch length optimization: %f\n", tr->likelihood);
 //printBranchLengths=TRUE;
 //Tree2String(tr->tree_string, tr, tr->start->back, printBranchLengths, FALSE, 0, 0, 0, SUMMARIZE_LH, 0,0);
 //fprintf(stderr, "%s\n", tr->tree_string);


 int foundNNI=TRUE;
 int nniRound=1;
 printf("Starting the first NNI search ... \n");
 getrusage(RUSAGE_SELF, &usage);
 start_tmp = usage.ru_utime;
 double curLH = tr->likelihood;
 while (TRUE) {
 double newLH = doNNISearch(tr);
 if (newLH == 0.0) {
 break;
 } else {
 printf("Likelihood after doing NNI round %d = %f\n", nniRound, newLH);
 //if (newLH - curLH < 0.001)
 //break;
 curLH = newLH;
 nniRound++;
 }

 }
 getrusage(RUSAGE_SELF, &usage);
 end_tmp = usage.ru_utime;
 printf("First NNI search finished in %f seconds \n", (double) end_tmp.tv_sec + end_tmp.tv_usec/1.0e6 - start_tmp.tv_sec - start_tmp.tv_usec/1.0e6);

 printf("Optimizing model parameters and branch lengths ... \n");
 getrusage(RUSAGE_SELF, &usage);
 start_tmp = usage.ru_utime;
 modOpt(tr, 0.1);
 getrusage(RUSAGE_SELF, &usage);
 end_tmp = usage.ru_utime;
 printf("Model parameters and branch lengths optimized in %f seconds \n", (double) end_tmp.tv_sec + end_tmp.tv_usec/1.0e6 - start_tmp.tv_sec - start_tmp.tv_usec/1.0e6);
 evaluateGeneric(tr, tr->start, FALSE);
 printf("Likelihood after model optimization: %f\n", tr->likelihood);
 Tree2String(tr->tree_string, tr, tr->start->back, printBranchLengths, TRUE, 0, 0, 0, SUMMARIZE_LH, 0,0);
 char resultFile[100];
 strcpy (resultFile, argv[1]);
 strcat(resultFile,".result");
 FILE *treefile = fopen(resultFile, "w");
 fprintf(treefile, "%s\n", tr->tree_string);
 getrusage(RUSAGE_SELF, &usage);
 end = usage.ru_utime;
 printf("CPU time = %f seconds \n", (double) end.tv_sec + end.tv_usec/1.0e6 - start.tv_sec - start.tv_usec/1.0e6);
 printf("Result tree written to %s \n", resultFile);
 //printTopology(tr, TRUE);
 fclose(treefile);
 return (0);
 }*/

void evalNNIForSubtree(pllInstance* tr, partitionList *pr, nodeptr p, NNIMOVE* nniList, int* numBran,
		int* numPosNNI, double curLH) {
	if (!isTip(p->number, tr->mxtips)) {
		int betterNNI = evalNNIForBran(tr, pr, p, nniList, numBran, curLH);
		if (betterNNI) {
			*numPosNNI = *numPosNNI + 1;
		}
		*numBran = *numBran + 2;
		nodeptr q = p->next;
		while (q != p) {
			evalNNIForSubtree(tr, pr, q->back, nniList, numBran, numPosNNI, curLH);
			q = q->next;
		}
	}
}
