/***************************************************************************
 *   Copyright (C) 2009 by BUI Quang Minh   *
 *   minh.bui@univie.ac.at   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef SUPERALIGNMENT_H
#define SUPERALIGNMENT_H

#include "alignment.h"


struct PartitionInfo {
	string name; // partition name
	string model_name; // model name
	string aln_file; // alignment file associated
	string sequence_type; // sequence type (DNA/AA/BIN)
	string position_spec; // position specification, e.g., "1-100\1 1-100\2"

	double cur_score; // current log-likelihood

	DoubleVector null_score; // log-likelihood of each branch collapsed to zero
	DoubleVector opt_score; // optimized log-likelihood for every branch
	DoubleVector nni1_score; // log-likelihood for 1st NNI for every branch
	DoubleVector nni2_score; // log-likelihood for 2nd NNI for every branch

	DoubleVector cur_brlen; // current branch lengths
	DoubleVector opt_brlen; // optimized branch lengths for every branch
	DoubleVector nni1_brlen; // branch length for 1st NNI for every branch
	DoubleVector nni2_brlen; // branch length for 2nd NNI for every branch

	double *mem_ptnlh; // total memory allocated for all pattern likelihood vectors
	double *cur_ptnlh; // current pattern likelihoods of the tree
	vector<double* > opt_ptnlh; // pattern likelihoods of tree where one branch is optimized
	vector<double* > nni1_ptnlh; // pattern likelihoods of 1st NNI tree
	vector<double* > nni2_ptnlh; // pattern likelihoods of 2nd NNI tree
};

class PhyloSuperTree;

/**
Super alignment representing presence/absence of sequences in
k partitions for a total of n sequences. It has the form:
		Site_1 Site_2 ... Site_k
Seq_1     1      0    ...   1
Seq_2     0      1    ...   0
...      ...
Seq_n     1      1    ...   0

Where (i,j)=1 means Seq_i is present in partition j, 0 otherwise

So data is binary.

	@author BUI Quang Minh <minh.bui@univie.ac.at>
*/

class SuperAlignment : public Alignment
{
public:
	/** constructor initialize from a supertree */
    SuperAlignment(PhyloSuperTree *super_tree);

	/** constructor initialize empty alignment */
    SuperAlignment();

    /** destructor */
    ~SuperAlignment();

    /** return that this is a super-alignment structure */
	virtual bool isSuperAlignment() { return true; }

	/**
	 * create taxa_index from super-alignment to sub-alignment
	 * @param part index of sub-alignment
	 */
	void linkSubAlignment(int part);

	/**
	 * @param pattern_index (OUT) vector of size = alignment length storing pattern index of all sites
	 * the index of sites in 2nd, 3rd,... genes have to be increased by the number of patterns in previous genes
	 * so that all indices are distinguishable
	*/
	virtual void getSitePatternIndex(IntVector &pattern_index);

	/**
	 * @param freq (OUT) vector of site-pattern frequencies for all sub-alignments
	*/
	virtual void getPatternFreq(IntVector &pattern_freq);

	/**
		Quit if some sequences contain only gaps or missing data
	*/
	//virtual void checkGappySeq(bool force_error = true);

	/**
		create a non-parametric bootstrap alignment by resampling sites within partitions
		@param aln input alignment
		@param pattern_freq (OUT) if not NULL, will store the resampled pattern frequencies
        @param spec bootstrap specification of the form "l1:b1,l2:b2,...,lk:bk"
            	to randomly draw b1 sites from the first l1 sites, etc. Note that l1+l2+...+lk
            	must equal m, where m is the alignment length. Otherwise, an error will occur.
            	If spec == NULL, a standard procedure is applied, i.e., randomly draw m sites.
	*/
	virtual void createBootstrapAlignment(Alignment *aln, IntVector* pattern_freq = NULL, const char *spec = NULL);

	/**
		resampling pattern frequency by a non-parametric bootstrap 
		@param pattern_freq (OUT) resampled pattern frequencies
        @param spec bootstrap specification, see above
	*/
	virtual void createBootstrapAlignment(IntVector &pattern_freq, const char *spec = NULL);

	/**
		resampling pattern frequency by a non-parametric bootstrap
		@param pattern_freq (OUT) resampled pattern frequencies
        @param spec bootstrap specification, see above
	*/
	virtual void createBootstrapAlignment(int *pattern_freq, const char *spec = NULL);

	/**
	 * shuffle alignment by randomizing the order of sites over all sub-alignments
	 */
	virtual void shuffleAlignment();

	/**
		compute the observed (Hamming) distance (number of different pairs of positions per site)
			between two sequences
		@param seq1 index of sequence 1
		@param seq2 index of sequence 2
		@return the observed distance between seq1 and seq2 (between 0.0 and 1.0)
	*/
	virtual double computeObsDist(int seq1, int seq2);

	/**
		compute the Juke-Cantor corrected distance between 2 sequences over all partitions
		@param seq1 index of sequence 1
		@param seq2 index of sequence 2		
		@return any distance between seq1 and seq2
	*/
	virtual double computeDist(int seq1, int seq2);

	/**
	 * print the super-alignment to a file
	 * @param filename
	 * @param append TRUE to append to this file, false to write new file
	 */
	void printCombinedAlignment(const char *filename, bool append = false);

	/**
	 * print all sub alignments into files with prefix, suffix is the charset name
	 * @param prefix prefix of output files
	 */
	void printSubAlignments(Params &params, vector<PartitionInfo> &part_info);

	/**
		@return unconstrained log-likelihood (without a tree)
	*/
	virtual double computeUnconstrainedLogL();

	/**
	 * @return proportion of missing data in super alignment
	 */
	double computeMissingData();

	/**
		actual partition alignments
	*/
	vector<Alignment*> partitions;

	/**
		matrix represents the index of taxon i in partition j, -1 if the taxon is not present
	*/
	vector<IntVector> taxa_index;

	/**
	 * concatenate subset of alignments
	 * @param ids IDs of sub-alignments
	 * @return concatenated alignment
	 */
    Alignment *concatenateAlignments(IntVector &ids);


};

#endif