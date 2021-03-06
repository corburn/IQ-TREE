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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iqtree_config.h>
#include "phylotree.h"
#include "phylosupertree.h"
#include "phylosupertreeplen.h"
#include "phyloanalysis.h"
#include "alignment.h"
#include "superalignment.h"
#include "iqtree.h"
#include "model/modelgtr.h"
#include "model/modeldna.h"
#include "myreader.h"
#include "model/rateheterogeneity.h"
#include "model/rategamma.h"
#include "model/rateinvar.h"
#include "model/rategammainvar.h"
//#include "modeltest_wrapper.h"
#include "model/modelprotein.h"
#include "model/modelbin.h"
#include "model/modelcodon.h"
#include "stoprule.h"

#include "mtreeset.h"
#include "mexttree.h"
#include "model/ratemeyerhaeseler.h"
#include "whtest_wrapper.h"
#include "model/partitionmodel.h"
#include "model/modelmixture.h"
//#include "guidedbootstrap.h"
#include "model/modelset.h"
#include "timeutil.h"
#include "upperbounds.h"


void reportReferences(Params &params, ofstream &out, string &original_model) {
	out << "To cite IQ-TREE please use:" << endl << endl
		<< "Lam-Tung Nguyen, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh (2015)" << endl
		<< "IQ-TREE: A fast and effective stochastic algorithm for estimating" << endl
		<< "maximum likelihood phylogenies. Mol. Biol. Evol., 32:268-274." << endl
        << "http://dx.doi.org/10.1093/molbev/msu300" << endl << endl;

	if (params.gbo_replicates)
	out << "Since you used ultrafast bootstrap (UFBoot) please also cite: " << endl << endl
		<< "Bui Quang Minh, Minh Anh Thi Nguyen, and Arndt von Haeseler (2013) Ultrafast" << endl
		<< "approximation for phylogenetic bootstrap. Mol. Biol. Evol., 30:1188-1195." << endl
        << "http://dx.doi.org/10.1093/molbev/mst024" << endl << endl;

    if (params.partition_file) 
    out << "Since you used partition models please also cite:" << endl << endl
        << "Olga Chernomor, Arndt von Haeseler, and Bui Quang Minh (2016) Terrace aware data" << endl
        << "structure for phylogenomic inference from supermatrices. Syst. Biol., in press." << endl
        << "http://dx.doi.org/10.1093/sysbio/syw037" << endl << endl;

}

void reportAlignment(ofstream &out, Alignment &alignment, int nremoved_seqs) {
	out << "Input data: " << alignment.getNSeq()+nremoved_seqs << " sequences with "
			<< alignment.getNSite() << " "
			<< ((alignment.seq_type == SEQ_BINARY) ?
					"binary" :
					((alignment.seq_type == SEQ_DNA) ? "nucleotide" :
					(alignment.seq_type == SEQ_PROTEIN) ? "amino-acid" :
					(alignment.seq_type == SEQ_CODON) ? "codon": "morphological"))
			<< " sites" 
            << endl << "Number of constant sites: "
			<< round(alignment.frac_const_sites * alignment.getNSite())
			<< " (= " << alignment.frac_const_sites * 100 << "% of all sites)"
            << endl << "Number of invariant (constant or ambiguous constant) sites: "
			<< round(alignment.frac_invariant_sites * alignment.getNSite())
			<< " (= " << alignment.frac_invariant_sites * 100 << "% of all sites)"
			<< endl << "Number of distinct site patterns: " << alignment.size() << endl
			<< endl;
}

void pruneModelInfo(vector<ModelInfo> &model_info, PhyloSuperTree *tree) {
	vector<ModelInfo> res_info;
	for (vector<PartitionInfo>::iterator it = tree->part_info.begin(); it != tree->part_info.end(); it++) {
		for (vector<ModelInfo>::iterator mit = model_info.begin(); mit != model_info.end(); mit++)
			if (mit->set_name == it->name)
				res_info.push_back(*mit);
	}
	model_info = res_info;

}

void reportModelSelection(ofstream &out, Params &params, vector<ModelInfo> &model_info, bool is_partitioned) {
	out << "Best-fit model according to "
		<< ((params.model_test_criterion == MTC_BIC) ? "BIC" :
			((params.model_test_criterion == MTC_AIC) ? "AIC" : "AICc")) << ": ";
	vector<ModelInfo>::iterator it;
	if (is_partitioned) {
		string set_name = "";
		for (it = model_info.begin(); it != model_info.end(); it++) {
			if (it->set_name != set_name) {
				if (set_name != "")
					out << ",";
				out << it->name << ":" << it->set_name;
				set_name = it->set_name;
			}
		}
	} else {
		out << model_info[0].name;
	}

	if (is_partitioned) {
		out << endl << endl << "List of best-fit models per partition:" << endl << endl;
	} else {
		out << endl << endl << "List of models sorted by "
			<< ((params.model_test_criterion == MTC_BIC) ? "BIC" :
				((params.model_test_criterion == MTC_AIC) ? "AIC" : "AICc"))
			<< " scores: " << endl << endl;
	}
	if (is_partitioned)
		out << "  ID  ";
	out << "Model             LogL          AIC      w-AIC      AICc     w-AICc       BIC      w-BIC" << endl;
	/*
	if (is_partitioned)
		out << "----------";

	out << "----------------------------------------------------------------------------------------" << endl;
	*/
	int setid = 1;
	for (it = model_info.begin(); it != model_info.end(); it++) {
		if (it->AIC_score == DBL_MAX) continue;
		if (it != model_info.begin() && it->set_name != (it-1)->set_name)
			setid++;
		if (is_partitioned && it != model_info.begin() && it->set_name == (it-1)->set_name)
			continue;
		if (is_partitioned) {
			out.width(4);
			out << right << setid << "  ";
		}
		out.width(15);
		out << left << it->name << " ";
		out.width(11);
		out << right << it->logl << " ";
		out.width(11);
		out	<< it->AIC_score << ((it->AIC_conf) ? " + " : " - ") << it->AIC_weight << " ";
		out.width(11);
		out << it->AICc_score << ((it->AICc_conf) ? " + " : " - ") << it->AICc_weight << " ";
		out.width(11);
		out << it->BIC_score  << ((it->BIC_conf) ? " + " : " - ") << it->BIC_weight;
		out << endl;
	}
	out << endl;
	out <<  "AIC, w-AIC   : Akaike information criterion scores and weights." << endl
		 << "AICc, w-AICc : Corrected AIC scores and weights." << endl
		 << "BIC, w-BIC   : Bayesian information criterion scores and weights." << endl << endl

		 << "Plus signs denote the 95% confidence sets." << endl
		 << "Minus signs denote significant exclusion." <<endl;
	out << endl;
}

void reportModel(ofstream &out, Alignment *aln, ModelSubst *m) {
	int i, j, k;
	assert(aln->num_states == m->num_states);
    double *rate_mat = new double[m->num_states * m->num_states];
    if (!m->isSiteSpecificModel())
        m->getRateMatrix(rate_mat);
    else
        ((ModelSet*)m)->front()->getRateMatrix(rate_mat);

	if (m->num_states <= 4) {
		out << "Rate parameter R:" << endl << endl;

		if (m->num_states > 4)
			out << fixed;
		if (m->isReversible()) {
			for (i = 0, k = 0; i < m->num_states - 1; i++)
				for (j = i + 1; j < m->num_states; j++, k++) {
					out << "  " << aln->convertStateBackStr(i) << "-" << aln->convertStateBackStr(j) << ": "
							<< rate_mat[k];
					if (m->num_states <= 4)
						out << endl;
					else if (k % 5 == 4)
						out << endl;
				}

		} else { // non-reversible model
			for (i = 0, k = 0; i < m->num_states; i++)
				for (j = 0; j < m->num_states; j++)
					if (i != j) {
						out << "  " << aln->convertStateBackStr(i) << "-" << aln->convertStateBackStr(j)
								<< ": " << rate_mat[k];
						if (m->num_states <= 4)
							out << endl;
						else if (k % 5 == 4)
							out << endl;
						k++;
					}

		}

		//if (tree.aln->num_states > 4)
		out << endl;
		out.unsetf(ios_base::fixed);
	} else if (aln->seq_type == SEQ_PROTEIN && m->getNDim() > 20) {
        assert(m->num_states == 20);
        out << "WARNING: This model has " << m->getNDim() + m->getNDimFreq() << " parameters that may be overfitting. Please use with caution!" << endl << endl;
        double full_mat[400];
        for (i = 0, k = 0; i < m->num_states - 1; i++)
            for (j = i + 1; j < m->num_states; j++, k++) {
                full_mat[i*m->num_states+j] = rate_mat[k];
            }
        out << "Substitution parameters (lower-diagonal) and state frequencies in PAML format (can be used as input for IQ-TREE): " << endl << endl;
        for (i = 1; i < m->num_states; i++) {
            for (j = 0; j < i; j++)
                out << "\t" << full_mat[j*m->num_states+i];
            out << endl;
        }
        double state_freq[20];
        m->getStateFrequency(state_freq);
        for (i = 0; i < m->num_states; i++)
            out << "\t" << state_freq[i];
        out << endl << endl;
    }
    
    delete[] rate_mat;

	out << "State frequencies: ";
	if (m->isSiteSpecificModel())
		out << "(site specific frequencies)" << endl << endl;
	else {
		if (!m->isReversible())
			out << "(inferred from Q matrix)" << endl;
		else
			switch (m->getFreqType()) {
			case FREQ_EMPIRICAL:
				out << "(empirical counts from alignment)" << endl;
				break;
			case FREQ_ESTIMATE:
				out << "(estimated with maximum likelihood)" << endl;
				break;
			case FREQ_USER_DEFINED:
				out << ((aln->seq_type == SEQ_PROTEIN) ? "(model)" : "(user-defined)") << endl;
				break;
			case FREQ_EQUAL:
				out << "(equal frequencies)" << endl;
				break;
			default:
				break;
			}
		out << endl;

		if (m->getFreqType() != FREQ_USER_DEFINED && m->getFreqType() != FREQ_EQUAL) {
			double *state_freqs = new double[m->num_states];
			m->getStateFrequency(state_freqs);
            int ncols=(aln->seq_type == SEQ_CODON) ? 4 : 1;
			for (i = 0; i < m->num_states; i++) {
				out << "  pi(" << aln->convertStateBackStr(i) << ") = " << state_freqs[i];
                if (i % ncols == ncols-1)
                    out << endl;
            }
			delete[] state_freqs;
			out << endl;
		}
		if (m->num_states <= 4) {
			// report Q matrix
			double *q_mat = new double[m->num_states * m->num_states];
			m->getQMatrix(q_mat);

			out << "Rate matrix Q:" << endl << endl;
			for (i = 0, k = 0; i < m->num_states; i++) {
				out << "  " << aln->convertStateBackStr(i);
				for (j = 0; j < m->num_states; j++, k++) {
					out << "  ";
					out.width(8);
					out << q_mat[k];
				}
				out << endl;
			}
			out << endl;
			delete[] q_mat;
		}
	}
}

void reportModel(ofstream &out, PhyloTree &tree) {
//	int i, j, k;
	int i;

	if (tree.getModel()->isMixture()) {
		out << "Mixture model of substitution: " << tree.getModelName() << endl;
//		out << "Full name: " << tree.getModelName() << endl;
		ModelMixture *mmodel = (ModelMixture*) tree.getModel();
		out << endl << "  No  Component      Rate    Weight   Parameters" << endl;
		i = 0;
		for (ModelMixture::iterator m = mmodel->begin(); m != mmodel->end(); m++, i++) {
			out.width(4);
			out << right << i+1 << "  ";
			out.width(12);
			out << left << (*m)->name << "  ";
			out.width(7);
			out << (*m)->total_num_subst << "  ";
			out.width(7);
			out << mmodel->prop[i] << "  " << (*m)->getNameParams() << endl;
//			out << "Model for mixture component "  << (m-mmodel->begin())+1 << ": " << (*m)->name << endl;
//			reportModel(out, tree.aln, *m);
		}
		out << endl;
	} else {
		out << "Model of substitution: " << tree.getModelName() << endl << endl;
		reportModel(out, tree.aln, tree.getModel());
	}
}

void reportRate(ofstream &out, PhyloTree &tree) {
	int i;
	RateHeterogeneity *rate_model = tree.getRate();
	out << "Model of rate heterogeneity: " << rate_model->full_name << endl;
	rate_model->writeInfo(out);

	if (rate_model->getNDiscreteRate() > 1 || rate_model->getPInvar() > 0.0) {
		out << endl << " Category  Relative_rate  Proportion" << endl;
		if (rate_model->getPInvar() > 0.0)
			out << "  0         0              " << rate_model->getPInvar()
					<< endl;
		int cats = rate_model->getNDiscreteRate();
		DoubleVector prop;
		if (rate_model->getGammaShape() > 0 || rate_model->getPtnCat(0) < 0) {
//			prop.resize(cats, (1.0 - rate_model->getPInvar()) / rate_model->getNRate());
			prop.resize(cats);
		for (i = 0; i < cats; i++)
			prop[i] = rate_model->getProp(i);
		} else {
			prop.resize(cats, 0.0);
			for (i = 0; i < tree.aln->getNPattern(); i++)
				prop[rate_model->getPtnCat(i)] += tree.aln->at(i).frequency;
			for (i = 0; i < cats; i++)
				prop[i] /= tree.aln->getNSite();
		}
		for (i = 0; i < cats; i++) {
			out << "  " << i + 1 << "         ";
			out.width(14);
			out << left << rate_model->getRate(i) << " " << prop[i];
			out << endl;
		}
		if (rate_model->isGammaRate()) {
			out << "Relative rates are computed as " << ((rate_model->isGammaRate() == GAMMA_CUT_MEDIAN) ? "MEDIAN" : "MEAN") <<
				" of the portion of the Gamma distribution falling in the category." << endl;
		}
	}
	/*
	 if (rate_model->getNDiscreteRate() > 1 || rate_model->isSiteSpecificRate())
	 out << endl << "See file " << rate_file << " for site-specific rates and categories" << endl;*/
	out << endl;
}

void reportTree(ofstream &out, Params &params, PhyloTree &tree, double tree_lh, double lh_variance, double main_tree) {
	double epsilon = 1.0 / tree.getAlnNSite();
	double totalLen = tree.treeLength();
	int df = tree.getModelFactory()->getNParameters();
	int ssize = tree.getAlnNSite();
	double AIC_score, AICc_score, BIC_score;
	computeInformationScores(tree_lh, df, ssize, AIC_score, AICc_score, BIC_score);
    
	out << "Log-likelihood of the tree: " << fixed << tree_lh;
    if (lh_variance > 0.0) 
        out << " (s.e. " << sqrt(lh_variance) << ")";
    out << endl;
    out	<< "Unconstrained log-likelihood (without tree): " << tree.aln->computeUnconstrainedLogL() << endl;

    out << "Number of free parameters (#branches + #model parameters): " << df << endl;
//    if (ssize > df) { 
//        if (ssize > 40*df)
//            out	<< "Akaike information criterion (AIC) score: " << AIC_score << endl;
//        else
//			out << "Corrected Akaike information criterion (AICc) score: " << AICc_score << endl;
//        
//		out << "Bayesian information criterion (BIC) score: " << BIC_score << endl;
//    } else 
    out	<< "Akaike information criterion (AIC) score: " << AIC_score << endl;
    out << "Corrected Akaike information criterion (AICc) score: " << AICc_score << endl;
    out << "Bayesian information criterion (BIC) score: " << BIC_score << endl;

    if (ssize <= df && main_tree) {
        
        out << endl
            << "**************************** WARNING ****************************" << endl
            << "Number of parameters (K, model parameters and branch lengths): " << df << endl
            << "Sample size (n, alignment length): " << ssize << endl << endl
            << "Given that K>=n, the parameter estimates might be inaccurate." << endl
            << "Thus, phylogenetic estimates should be interpreted with caution." << endl << endl 

            << "Ideally, it is desirable that n >> K. When selecting optimal models," << endl
            << "1. use AIC or BIC if n > 40K;" << endl 
            << "2. use AICc or BIC if 40K >= n > K;" << endl 
            << "3. be extremely cautious if n <= K" << endl << endl

            << "To improve the situation (3), consider the following options:" << endl
            << "  1. Increase the sample size (n)" << endl
            << "  2. Decrease the number of parameters (K) to be estimated. If" << endl
            << "     possible:" << endl
            << "     a. Remove the least important sequences from the alignment" << endl
            << "     b. Specify some of the parameter values for the substitution"<< endl 
            << "        model (e.g., the nucleotide or amino acid frequencies)" << endl
            << "     c. Specify some of the parameter values for the rates-across-" << endl
            << "        sites model (e.g., the shape parameter for the discrete" << endl
            << "        Gamma distribution, the proportion of invariable sites, or" << endl
            << "        the rates of change for different rate categories under" << endl
            << "        the FreeRate model)" << endl << endl
            << "Reference:" << endl
            << "Burnham KR, Anderson DR (2002). Model Selection and Multimodel" << endl
            << "Inference: A Practical Information-Theoretic Approach. Springer," << endl
            << "New York." << endl 
            << "************************ END OF WARNING ***********************" << endl;
    }
    out << endl;
    
	out << "Total tree length (sum of branch lengths): " << totalLen << endl;
	double totalLenInternal = tree.treeLengthInternal(epsilon);
	out << "Sum of internal branch lengths: " << totalLenInternal << " (" << totalLenInternal*100.0 / totalLen << "% of tree length)" << endl;
//	out << "Sum of internal branch lengths divided by total tree length: "
//			<< totalLenInternal / totalLen << endl;
	out << endl;
	//out << "ZERO BRANCH EPSILON = " << epsilon << endl;
	int zero_internal_branches = tree.countZeroInternalBranches(NULL, NULL, epsilon);
	if (zero_internal_branches > 0) {
		//int zero_internal_branches = tree.countZeroInternalBranches(NULL, NULL, epsilon);
		/*
		out << "WARNING: " << zero_branches
				<< " branches of near-zero lengths (<" << epsilon << ") and should be treated with caution!"
				<< endl;
		*/
		out << "WARNING: " << zero_internal_branches
				<< " near-zero internal branches (<" << epsilon << ") should be treated with caution"
				<< endl;
		/*
		cout << endl << "WARNING: " << zero_branches
				<< " branches of near-zero lengths (<" << epsilon << ") and should be treated with caution!"
				<< endl;
		*/
		out << "         Such branches are denoted by '**' in the figure below"
				<< endl << endl;
	}
	int long_branches = tree.countLongBranches(NULL, NULL, params.max_branch_length-0.2);
	if (long_branches > 0) {
		//stringstream sstr;
		out << "WARNING: " << long_branches << " too long branches (>" 
            << params.max_branch_length-0.2 << ") should be treated with caution!" << endl;
		//out << sstr.str();
		//cout << sstr.str();
	}

			//<< "Total tree length: " << tree.treeLength() << endl << endl
	tree.sortTaxa();
    out << "NOTE: Tree is UNROOTED although outgroup taxon '" << tree.root->name << "' is drawn at root" << endl;

    if (tree.isSuperTree() && params.partition_type == 0)
        out	<< "NOTE: Branch lengths are weighted average over all partitions" << endl
            << "      (weighted by the number of sites in the partitions)" << endl;

    bool is_codon = tree.aln->seq_type == SEQ_CODON;
    if (tree.isSuperTree()) {
        PhyloSuperTree *stree = (PhyloSuperTree*) &tree;
        is_codon = true;
        for (PhyloSuperTree::iterator sit = stree->begin(); sit != stree->end(); sit++)
            if ((*sit)->aln->seq_type != SEQ_CODON) {
                is_codon = false;
                break;
            }
    }
    if (is_codon)
		out << endl << "NOTE: Branch lengths are interpreted as number of nucleotide substitutions per codon site!"
				<< endl << "      Rescale them by 1/3 if you want to have #nt substitutions per nt site" << endl;
    if (main_tree) 
    if (params.aLRT_replicates > 0 || params.gbo_replicates || (params.num_bootstrap_samples && params.compute_ml_tree)) {
        out << "Numbers in parentheses are ";
        if (params.aLRT_replicates > 0) {
            out << "SH-aLRT support (%)";
            if (params.localbp_replicates)
                out << " / local bootstrap support (%)";
        }
        if (params.aLRT_test)
            out << " / parametric aLRT support";
        if (params.aBayes_test)
            out << " / aBayes support";
        if (params.num_bootstrap_samples && params.compute_ml_tree) {
            if (params.aLRT_replicates > 0 || params.aLRT_test || params.aBayes_test)
                out << " /";
            out << " standard bootstrap support (%)";
        }
        if (params.gbo_replicates) {
            if (params.aLRT_replicates > 0 || params.aLRT_test || params.aBayes_test)
                out << " /";
            out << " ultrafast bootstrap support (%)";
        }
        out << endl;
    }
    out << endl;

	//tree.setExtendedFigChar();
	tree.drawTree(out, WT_BR_SCALE, epsilon);
        
    out << "Tree in newick format:" << endl << endl;

	tree.printTree(out, WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA);

	out << endl << endl;
}

void reportCredits(ofstream &out) {
	out << "CREDITS" << endl << "-------" << endl << endl
			<< "Some parts of the code were taken from the following packages/libraries:"
			<< endl << endl
			<< "Schmidt HA, Strimmer K, Vingron M, and von Haeseler A (2002)" << endl
			<< "TREE-PUZZLE: maximum likelihood phylogenetic analysis using quartets" << endl
			<< "and parallel computing. Bioinformatics, 18(3):502-504." << endl << endl

			//<< "The source code to construct the BIONJ tree were taken from BIONJ software:"
			//<< endl << endl
			<< "Gascuel O (1997) BIONJ: an improved version of the NJ algorithm" << endl
			<< "based on a simple model of sequence data. Mol. Bio. Evol., 14:685-695." << endl << endl

			//<< "The Nexus file parser was taken from the Nexus Class Library:"
			//<< endl << endl
			<< "Paul O. Lewis (2003) NCL: a C++ class library for interpreting data files in" << endl
			<< "NEXUS format. Bioinformatics, 19(17):2330-2331." << endl << endl

			<< "Mascagni M and Srinivasan A (2000) Algorithm 806: SPRNG: A Scalable Library" << endl
			<< "for Pseudorandom Number Generation. ACM Transactions on Mathematical Software," << endl
			<< "26: 436-461." << endl << endl

			<< "Guennebaud G, Jacob B, et al. (2010) Eigen v3. http://eigen.tuxfamily.org" << endl << endl;
			/*
			<< "The Modeltest 3.7 source codes were taken from:" << endl << endl
			<< "David Posada and Keith A. Crandall (1998) MODELTEST: testing the model of"
			<< endl << "DNA substitution. Bioinformatics, 14(9):817-8." << endl
			*/
}

/***********************************************************
 * CREATE REPORT FILE
 ***********************************************************/
extern StringIntMap pllTreeCounter;

void exhaustiveSearchGAMMAInvar(Params &params, IQTree &iqtree);

void searchGAMMAInvarByRestarting(IQTree &iqtree);

void computeLoglFromUserInputGAMMAInvar(Params &params, IQTree &iqtree);

void printOutfilesInfo(Params &params, string &original_model, IQTree &tree) {

	cout << endl << "Analysis results written to: " << endl;
    if (!(params.suppress_output_flags & OUT_IQTREE))
		cout<< "  IQ-TREE report:                " << params.out_prefix << ".iqtree"
			<< endl;
	if (params.compute_ml_tree) {
        if (!(params.suppress_output_flags & OUT_TREEFILE)) {
            if (original_model.find("ONLY") == string::npos)
                cout << "  Maximum-likelihood tree:       " << params.out_prefix << ".treefile" << endl;
            else
                cout << "  Tree used for model selection: " << params.out_prefix << ".treefile" << endl;
        }
		if (params.snni && params.write_local_optimal_trees) {
			cout << "  Locally optimal trees (" << tree.candidateTrees.getNumLocalOptTrees() << "):    " << params.out_prefix << ".suboptimal_trees" << endl;
		}
	}
	if (!params.user_file && params.start_tree == STT_BIONJ) {
		cout << "  BIONJ tree:                    " << params.out_prefix << ".bionj"
				<< endl;
	}
	if (!params.dist_file) {
		//cout << "  Juke-Cantor distances:    " << params.out_prefix << ".jcdist" << endl;
		if (params.compute_ml_dist)
		cout << "  Likelihood distances:          " << params.out_prefix
					<< ".mldist" << endl;
		if (params.print_conaln)
		cout << "  Concatenated alignment:        " << params.out_prefix
					<< ".conaln" << endl;
	}
	if (original_model.find("TEST") != string::npos && tree.isSuperTree()) {
		cout << "  Best partitioning scheme:      " << params.out_prefix << ".best_scheme.nex" << endl;
		bool raxml_format_printed = true;

		for (vector<PartitionInfo>::iterator it = ((PhyloSuperTree*)&tree)->part_info.begin();
				it != ((PhyloSuperTree*)&tree)->part_info.end(); it++)
			if (!it->aln_file.empty()) {
				raxml_format_printed = false;
				break;
			}
		if (raxml_format_printed)
			 cout << "           in RAxML format:      " << params.out_prefix << ".best_scheme" << endl;
	}
	if (tree.getRate()->getGammaShape() > 0 && params.print_site_rate)
		cout << "  Gamma-distributed rates:       " << params.out_prefix << ".rate"
				<< endl;

	if ((tree.getRate()->isSiteSpecificRate() || tree.getRate()->getPtnCat(0) >= 0) && params.print_site_rate)
		cout << "  Site-rates by MH model:        " << params.out_prefix << ".rate"
				<< endl;

	if (params.print_site_lh)
		cout << "  Site log-likelihoods:          " << params.out_prefix << ".sitelh"
				<< endl;

	if (params.print_partition_lh)
		cout << "  Partition log-likelihoods:     " << params.out_prefix << ".partlh"
				<< endl;

	if (params.print_site_prob)
		cout << "  Site probability per rate/mix: " << params.out_prefix << ".siteprob"
				<< endl;

	if (params.write_intermediate_trees)
		cout << "  All intermediate trees:        " << params.out_prefix << ".treels"
				<< endl;

	if (params.gbo_replicates) {
		cout << endl << "Ultrafast bootstrap approximation results written to:" << endl
			 << "  Split support values:          " << params.out_prefix << ".splits.nex" << endl
			 << "  Consensus tree:                " << params.out_prefix << ".contree" << endl;
		if (params.print_ufboot_trees)
		cout << "  UFBoot trees:                  " << params.out_prefix << ".ufboot" << endl;

	}

	if (params.treeset_file) {
		cout << "  Evaluated user trees:          " << params.out_prefix << ".trees" << endl;

		if (params.print_tree_lh) {
		cout << "  Tree log-likelihoods:          " << params.out_prefix << ".treelh" << endl;
		}
	}
    	if (params.lmap_num_quartets >= 0) {
		cout << "  Likelihood mapping plot (SVG): " << params.out_prefix << ".lmap.svg" << endl;
		cout << "  Likelihood mapping plot (EPS): " << params.out_prefix << ".lmap.eps" << endl;
	}
    if (!(params.suppress_output_flags & OUT_LOG))
    	cout << "  Screen log file:               " << params.out_prefix << ".log" << endl;
	/*	if (original_model == "WHTEST")
	 cout <<"  WH-TEST report:           " << params.out_prefix << ".whtest" << endl;*/
	cout << endl;

}


void reportPhyloAnalysis(Params &params, string &original_model,
		IQTree &tree, vector<ModelInfo> &model_info) {
        
    if (params.suppress_output_flags & OUT_IQTREE) {
        printOutfilesInfo(params, original_model, tree);
        return;
    }
        
	if (params.count_trees) {
		// addon: print #distinct trees
		cout << endl << "NOTE: " << pllTreeCounter.size() << " distinct trees evaluated during whole tree search" << endl;

		IntVector counts;
		for (StringIntMap::iterator i = pllTreeCounter.begin(); i != pllTreeCounter.end(); i++) {
			if (i->second > counts.size())
				counts.resize(i->second+1, 0);
			counts[i->second]++;
		}
		for (IntVector::iterator i2 = counts.begin(); i2 != counts.end(); i2++) {
		    if (*i2 != 0) {
				cout << "#Trees occurring " << (i2-counts.begin()) << " times: " << *i2 << endl;
		    }
		}
	}
	string outfile = params.out_prefix;

	outfile += ".iqtree";
	try {
		ofstream out;
		out.exceptions(ios::failbit | ios::badbit);
		out.open(outfile.c_str());
		out << "IQ-TREE " << iqtree_VERSION_MAJOR << "." << iqtree_VERSION_MINOR
				<< "." << iqtree_VERSION_PATCH << " built " << __DATE__ << endl
				<< endl;
		if (params.partition_file)
			out << "Partition file name: " << params.partition_file << endl;
		if (params.aln_file)
			out << "Input file name: " << params.aln_file << endl;

		if (params.user_file)
			out << "User tree file name: " << params.user_file << endl;
		out << "Type of analysis: ";
        if (original_model.find("TEST") != string::npos && original_model.find("ONLY") != string::npos) {
            out << "model selection";
        } else {
            if (params.compute_ml_tree)
                out << "tree reconstruction";
            if (params.num_bootstrap_samples > 0) {
                if (params.compute_ml_tree)
                    out << " + ";
                out << "non-parametric bootstrap (" << params.num_bootstrap_samples
                        << " replicates)";
            }
            if (params.gbo_replicates > 0) {
                out << " + ultrafast bootstrap (" << params.gbo_replicates << " replicates)";
            }
        }
		out << endl;
		out << "Random seed number: " << params.ran_seed << endl << endl;
		out << "REFERENCES" << endl << "----------" << endl << endl;
		reportReferences(params, out, original_model);

		out << "SEQUENCE ALIGNMENT" << endl << "------------------" << endl
				<< endl;
		if (tree.isSuperTree()) {
			out << "Input data: " << tree.aln->getNSeq()+tree.removed_seqs.size() << " taxa with "
					<< tree.aln->getNSite() << " partitions and "
					<< tree.getAlnNSite() << " total sites ("
					<< ((SuperAlignment*)tree.aln)->computeMissingData()*100 << "% missing data)" << endl << endl;

			PhyloSuperTree *stree = (PhyloSuperTree*) &tree;
			int namelen = stree->getMaxPartNameLength();
			int part;
			out.width(max(namelen+6,10));
			out << left << "  ID  Name" << "  Type  #Seqs  #Sites  #Patterns  #Const_Sites" << endl;
			//out << string(namelen+54, '-') << endl;
			part = 0;
			for (PhyloSuperTree::iterator it = stree->begin(); it != stree->end(); it++, part++) {
				//out << "FOR PARTITION " << stree->part_info[part].name << ":" << endl << endl;
				//reportAlignment(out, *((*it)->aln));
				out.width(4);
				out << right << part+1 << "  ";
				out.width(max(namelen,4));
				out << left << stree->part_info[part].name << "  ";
				out.width(6);
				switch ((*it)->aln->seq_type) {
				case SEQ_BINARY: out << "BIN"; break;
				case SEQ_CODON: out << "CODON"; break;
				case SEQ_DNA: out << "DNA"; break;
				case SEQ_MORPH: out << "MORPH"; break;
				case SEQ_MULTISTATE: out << "TINA"; break;
				case SEQ_PROTEIN: out << "AA"; break;
				case SEQ_UNKNOWN: out << "???"; break;
				}
				out.width(5);
				out << right << (*it)->aln->getNSeq() << "  ";
				out.width(6);
				out << (*it)->aln->getNSite() << "  ";
				out.width(6);
				out << (*it)->aln->getNPattern() << "      ";
				out << round((*it)->aln->frac_const_sites*100) << "%" << endl;
			}
			out << endl;
		} else
			reportAlignment(out, *(tree.aln), tree.removed_seqs.size());

		out.precision(4);
		out << fixed;

		if (!model_info.empty()) {
			out << "MODEL SELECTION" << endl << "---------------" << endl << endl;
			if (tree.isSuperTree())
				pruneModelInfo(model_info, (PhyloSuperTree*)&tree);
			reportModelSelection(out, params, model_info, tree.isSuperTree());
		}

		out << "SUBSTITUTION PROCESS" << endl << "--------------------" << endl
				<< endl;
		if (tree.isSuperTree()) {
			if(params.partition_type)
				out	<< "Proportional partition model with joint branch lengths and separate models between partitions" << endl << endl;
			else
				out	<< "Full partition model with separate branch lengths and models between partitions" << endl << endl;
			PhyloSuperTree *stree = (PhyloSuperTree*) &tree;
			PhyloSuperTree::iterator it;
			int part;
			if(params.partition_type)
				out << "  ID  Model           Speed  Parameters" << endl;
			else
				out << "  ID  Model         TreeLen  Parameters" << endl;
			//out << "-------------------------------------" << endl;
			for (it = stree->begin(), part = 0; it != stree->end(); it++, part++) {
				out.width(4);
				out << right << (part+1) << "  ";
				out.width(14);
				if(params.partition_type)
					out << left << (*it)->getModelName() << " " << stree->part_info[part].part_rate  << "  " << (*it)->getModelNameParams() << endl;
				else
					out << left << (*it)->getModelName() << " " << (*it)->treeLength() << "  " << (*it)->getModelNameParams() << endl;
			}
			out << endl;
			/*
			for (it = stree->begin(), part = 0; it != stree->end(); it++, part++) {
				reportModel(out, *(*it));
				reportRate(out, *(*it));
			}*/
		} else {
			reportModel(out, tree);
			reportRate(out, tree);
		}

    		if (params.lmap_num_quartets >= 0) {
			tree.reportLikelihoodMapping(out);
		}


		/*
		out << "RATE HETEROGENEITY" << endl << "------------------" << endl
				<< endl;
		if (tree.isSuperTree()) {
			PhyloSuperTree *stree = (PhyloSuperTree*) &tree;
			int part = 0;
			for (PhyloSuperTree::iterator it = stree->begin();
					it != stree->end(); it++, part++) {
				out << "FOR PARTITION " << stree->part_info[part].name << ":"
						<< endl << endl;
				reportRate(out, *(*it));
			}
		} else
			reportRate(out, tree);
		*/
		// Bootstrap analysis:
		//Display as outgroup: a

		if (original_model == "WHTEST") {
			out << "TEST OF MODEL HOMOGENEITY" << endl
					<< "-------------------------" << endl << endl;
			out << "Delta of input data:                 "
					<< params.whtest_delta << endl;
			out << ".95 quantile of Delta distribution:  "
					<< params.whtest_delta_quantile << endl;
			out << "Number of simulations performed:     "
					<< params.whtest_simulations << endl;
			out << "P-value:                             "
					<< params.whtest_p_value << endl;
			if (params.whtest_p_value < 0.05) {
				out
						<< "RESULT: Homogeneity assumption is rejected (p-value cutoff 0.05)"
						<< endl;
			} else {
				out
						<< "RESULT: Homogeneity assumption is NOT rejected (p-value cutoff 0.05)"
						<< endl;
			}
			out << endl << "*** For this result please cite:" << endl << endl;
			out
					<< "G. Weiss and A. von Haeseler (2003) Testing substitution models"
					<< endl
					<< "within a phylogenetic tree. Mol. Biol. Evol, 20(4):572-578"
					<< endl << endl;
		}
/*
		out << "TREE SEARCH" << endl << "-----------" << endl << endl
				<< "Stopping rule: "
				<< ((params.stop_condition == SC_STOP_PREDICT) ? "Yes" : "No")
				<< endl << "Number of iterations: "
				<< tree.stop_rule.getNumIterations() << endl
				<< "Probability of deleting sequences: " << params.p_delete
				<< endl << "Number of representative leaves: "
				<< params.k_representative << endl
				<< "NNI log-likelihood cutoff: " << tree.getNNICutoff() << endl
				<< endl;
*/
		if (params.compute_ml_tree) {
			if (original_model.find("ONLY") != string::npos) {
				out << "TREE USED FOR MODEL SELECTION" << endl
					<< "-----------------------------" << endl << endl;
            } else if (params.min_iterations == 0) {
                if (params.user_file)
                    out << "USER TREE" << endl
                        << "---------" << endl << endl;
                else
                    out << "STARTING TREE" << endl
                        << "-------------" << endl << endl;                
            } else { 
				out << "MAXIMUM LIKELIHOOD TREE" << endl
					<< "-----------------------" << endl << endl;
            }

			tree.setRootNode(params.root);
            
            if (params.gbo_replicates) {
                if (tree.boot_consense_logl > tree.candidateTrees.getBestScore() + 0.1) {
                    out << endl << "**NOTE**: Consensus tree has higher likelihood than ML tree found! Please use consensus tree below." << endl;
                }
            }

			reportTree(out, params, tree, tree.candidateTrees.getBestScore(), tree.logl_variance, true);

			if (tree.isSuperTree() && verbose_mode >= VB_MED) {
				PhyloSuperTree *stree = (PhyloSuperTree*) &tree;
//				stree->mapTrees();
//				int empty_branches = stree->countEmptyBranches();
//				if (empty_branches) {
//					stringstream ss;
//					ss << empty_branches << " branches in the overall tree with no phylogenetic information due to missing data!";
//					outWarning(ss.str());
//				}
				
				int part = 0;
				for (PhyloSuperTree::iterator it = stree->begin();
						it != stree->end(); it++, part++) {
					out << "FOR PARTITION " << stree->part_info[part].name
							<< ":" << endl << endl;
					string root_name;
					if (params.root)
						root_name = params.root;
					else
						root_name = (*it)->aln->getSeqName(0);
					(*it)->root = (*it)->findNodeName(root_name);
					assert((*it)->root);
//					reportTree(out, params, *(*it), (*it)->computeLikelihood(), (*it)->computeLogLVariance(), false);
					reportTree(out, params, *(*it), stree->part_info[part].cur_score, 0.0, false);
				}
			}

		}
		/*
		 if (params.write_intermediate_trees) {
		 out << endl << "CONSENSUS OF INTERMEDIATE TREES" << endl << "-----------------------" << endl << endl
		 << "Number of intermediate trees: " << tree.stop_rule.getNumIterations() << endl
		 << "Split threshold: " << params.split_threshold << endl
		 << "Burn-in: " << params.tree_burnin << endl << endl;
		 }*/

		if (params.consensus_type == CT_CONSENSUS_TREE) {
			out << "CONSENSUS TREE" << endl << "--------------" << endl << endl;
			out << "Consensus tree is constructed from "
					<< (params.num_bootstrap_samples ? params.num_bootstrap_samples : params.gbo_replicates)
					<< " bootstrap trees";
            if (params.gbo_replicates) {
                out << endl << "Log-likelihood of consensus tree: " << tree.boot_consense_logl;
            }
			string con_file = params.out_prefix;
			con_file += ".contree";

            IntVector rfdist;
            tree.computeRFDist(con_file.c_str(), rfdist);
            out << endl << "Robinson-Foulds distance between ML tree and consensus tree: " << rfdist[0] << endl;
            
            out << endl << "Branches with bootstrap support >"
					<< floor(params.split_threshold * 1000) / 10 << "% are kept";
			if (params.split_threshold == 0.0)
				out << " (extended consensus)";
			if (params.split_threshold == 0.5)
				out << " (majority-rule consensus)";
			if (params.split_threshold >= 0.99)
				out << " (strict consensus)";

			out << endl << "Branch lengths are optimized by maximum likelihood on original alignment" << endl;
			out << "Numbers in parentheses are bootstrap supports (%)" << endl << endl;

			bool rooted = false;
			MTree contree;
			contree.readTree(con_file.c_str(), rooted);
			contree.drawTree(out, WT_BR_SCALE);
			out << endl << "Consensus tree in newick format: " << endl << endl;
			contree.printTree(out);
			out << endl << endl;
//			tree.freeNode();
//			tree.root = NULL;
//			tree.readTree(con_file.c_str(), rooted);
//			if (removed_seqs.size() > 0) {
//				tree.reinsertIdenticalSeqs(tree.aln, removed_seqs, twin_seqs);
//			}
//			tree.setAlignment(tree.aln);

			// bug fix
//			if ((tree.sse == LK_EIGEN || tree.sse == LK_EIGEN_SSE) && !tree.isBifurcating()) {
//				cout << "NOTE: Changing to old kernel as consensus tree is multifurcating" << endl;
//				tree.changeLikelihoodKernel(LK_SSE);
//			}

//			tree.initializeAllPartialLh();
//			tree.fixNegativeBranch(false);
//			if (tree.isSuperTree())
//				((PhyloSuperTree*) &tree)->mapTrees();
//			tree.optimizeAllBranches();
//			tree.printTree(con_file.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA);
//			tree.sortTaxa();
//			tree.drawTree(out, WT_BR_SCALE);
//			out << endl << "Consensus tree in newick format: " << endl << endl;
//			tree.printResultTree(out);
//			out << endl << endl;
		}


		/* evaluate user trees */
		vector<TreeInfo> info;
		IntVector distinct_trees;
		if (params.treeset_file) {
			evaluateTrees(params, &tree, info, distinct_trees);
			out.precision(4);
            out.setf(ios_base::fixed);

			out << endl << "USER TREES" << endl << "----------" << endl << endl;
			out << "See " << params.out_prefix << ".trees for trees with branch lengths." << endl << endl;
			if (params.topotest_replicates && info.size() > 1) {
                out << "Tree      logL    deltaL  bp-RELL    p-KH     p-SH    ";
				if (params.do_weighted_test)
					out << "p-WKH    p-WSH    ";
                out << "c-ELW";
                if (params.do_au_test) 
                    out << "     p-AU";
                    
                out << endl << "------------------------------------------------------------------";
                if (params.do_weighted_test) 
					out << "------------------";
                if (params.do_au_test)
                    out << "-------";
                out << endl;
			} else {
				out << "Tree      logL    deltaL" << endl;
				out << "-------------------------" << endl;

			}
			double maxL = -DBL_MAX;
			int tid, orig_id;
			for (tid = 0; tid < info.size(); tid++)
				if (info[tid].logl > maxL) maxL = info[tid].logl;
			for (orig_id = 0, tid = 0; orig_id < distinct_trees.size(); orig_id++) {
				out.width(3);
				out << right << orig_id+1 << " ";
				if (distinct_trees[orig_id] >= 0) {
					out << " = tree " << distinct_trees[orig_id]+1 << endl;
					continue;
				}
				out.precision(3);
				out.width(12);
				out << info[tid].logl << " ";
				out.width(7);
				out << maxL - info[tid].logl;
				if (!params.topotest_replicates || info.size() <= 1) {
					out << endl;
					tid++;
					continue;
				}
				out.precision(4);
				out << "  ";
				out.width(6);
				out << info[tid].rell_bp;
				if (info[tid].rell_confident)
					out << " + ";
				else
					out << " - ";
				out.width(6);
				out << right << info[tid].kh_pvalue;
				if (info[tid].kh_pvalue < 0.05)
					out << " - ";
				else
					out << " + ";
				out.width(6);
				out << right << info[tid].sh_pvalue;
				if (info[tid].sh_pvalue < 0.05)
					out << " - ";
				else
					out << " + ";
                
				if (params.do_weighted_test) {
					out.width(6);
					out << right << info[tid].wkh_pvalue;
					if (info[tid].wkh_pvalue < 0.05)
						out << " - ";
					else
						out << " + ";
					out.width(6);
					out << right << info[tid].wsh_pvalue;
					if (info[tid].wsh_pvalue < 0.05)
						out << " - ";
					else
						out << " + ";
				}
				out.width(6);
				out << info[tid].elw_value;
				if (info[tid].elw_confident)
					out << " + ";
				else
					out << " - ";

                if (params.do_au_test) {
                    out.width(6);
                    out << right << info[tid].au_pvalue;
                    if (info[tid].au_pvalue < 0.05)
                        out << " - ";
                    else
                        out << " + ";
                }

				out << endl;
				tid++;
			}
			out << endl;

			if (params.topotest_replicates) {
				out <<  "deltaL  : logL difference from the maximal logl in the set." << endl
					 << "bp-RELL : bootstrap proportion using RELL method (Kishino et al. 1990)." << endl
					 << "p-KH    : p-value of one sided Kishino-Hasegawa test (1989)." << endl
					 << "p-SH    : p-value of Shimodaira-Hasegawa test (2000)." << endl;
				if (params.do_weighted_test) {
					out << "p-WKH   : p-value of weighted KH test." << endl
					 << "p-WSH   : p-value of weighted SH test." << endl;
				}
				out	 << "c-ELW   : Expected Likelihood Weight (Strimmer & Rambaut 2002)." << endl;
                if (params.do_au_test) {
					out << "p-AU    : p-value of approximately unbiased (AU) test (Shimodaira, 2002)." << endl;
                }
                out  << endl
					 << "Plus signs denote the 95% confidence sets." << endl
					 << "Minus signs denote significant exclusion."  << endl
					 << "All tests performed "
					 << params.topotest_replicates << " resamplings using the RELL method."<<endl;
			}
			out << endl;
		}


		time_t cur_time;
		time(&cur_time);

		char *date_str;
		date_str = ctime(&cur_time);
		out.unsetf(ios_base::fixed);
		out << "TIME STAMP" << endl << "----------" << endl << endl
				<< "Date and time: " << date_str << "Total CPU time used: "
				<< (double) params.run_time << " seconds (" << convert_time(params.run_time) << ")" << endl
				<< "Total wall-clock time used: " << getRealTime() - params.start_real_time
				<< " seconds (" << convert_time(getRealTime() - params.start_real_time) << ")" << endl << endl;

		//reportCredits(out); // not needed, now in the manual
		out.close();

	} catch (ios::failure) {
		outError(ERR_WRITE_OUTPUT, outfile);
	}
    
    printOutfilesInfo(params, original_model, tree);
}

void checkZeroDist(Alignment *aln, double *dist) {
	int ntaxa = aln->getNSeq();
	IntVector checked;
	checked.resize(ntaxa, 0);
	int i, j;
	for (i = 0; i < ntaxa - 1; i++) {
		if (checked[i])
			continue;
		string str = "";
		bool first = true;
		for (j = i + 1; j < ntaxa; j++)
			if (dist[i * ntaxa + j] <= Params::getInstance().min_branch_length) {
				if (first)
					str = "ZERO distance between sequences "
							+ aln->getSeqName(i);
				str += ", " + aln->getSeqName(j);
				checked[j] = 1;
				first = false;
			}
		checked[i] = 1;
		if (str != "")
			outWarning(str);
	}
}


void printAnalysisInfo(int model_df, IQTree& iqtree, Params& params) {
//	if (!params.raxmllib) {
	cout << "Model of evolution: ";
	if (iqtree.isSuperTree()) {
		cout << iqtree.getModelName() << " (" << model_df << " free parameters)" << endl;
	} else {
		cout << iqtree.getModelName() << " with ";
		switch (iqtree.getModel()->getFreqType()) {
		case FREQ_EQUAL:
			cout << "equal";
			break;
		case FREQ_EMPIRICAL:
			cout << "counted";
			break;
		case FREQ_USER_DEFINED:
			cout << "user-defined";
			break;
		case FREQ_ESTIMATE:
			cout << "optimized";
			break;
		case FREQ_CODON_1x4:
			cout << "counted 1x4";
			break;
		case FREQ_CODON_3x4:
			cout << "counted 3x4";
			break;
		case FREQ_CODON_3x4C:
			cout << "counted 3x4-corrected";
			break;
		default:
			outError("Wrong specified state frequencies");
		}
		cout << " frequencies (" << model_df << " free parameters)" << endl;
	}
	cout << "Fixed branch lengths: "
			<< ((params.fixed_branch_length) ? "Yes" : "No") << endl;

	if (params.min_iterations > 0) {
	    cout << "Tree search algorithm: " << (params.snni ? "Stochastic nearest neighbor interchange" : "IQPNNI") << endl;
	    cout << "Termination condition: ";
	    if (params.stop_condition == SC_REAL_TIME) {
	        cout << "after " << params.maxtime << " minutes" << endl;
	    } else if (params.stop_condition == SC_UNSUCCESS_ITERATION) {
	        cout << "after " << params.unsuccess_iteration << " unsuccessful iterations" << endl;
	    } else if (params.stop_condition == SC_FIXED_ITERATION) {
	            cout << params.min_iterations << " iterations" << endl;
	    } else if(params.stop_condition == SC_WEIBULL) {
	            cout << "predicted in [" << params.min_iterations << ","
	                    << params.max_iterations << "] (confidence "
	                    << params.stop_confidence << ")" << endl;
	    } else if (params.stop_condition == SC_BOOTSTRAP_CORRELATION) {
	    	cout << "min " << params.min_correlation << " correlation coefficient" << endl;
	    }

	    if (!params.snni) {
	        cout << "Number of representative leaves  : " << params.k_representative << endl;
	        cout << "Probability of deleting sequences: " << iqtree.getProbDelete() << endl;
	        cout << "Number of leaves to be deleted   : " << iqtree.getDelete() << endl;
	        cout << "Important quartets assessed on: "
	                << ((params.iqp_assess_quartet == IQP_DISTANCE) ?
	                        "Distance" : ((params.iqp_assess_quartet == IQP_PARSIMONY) ? "Parsimony" : "Bootstrap"))
	                << endl;
	    }
	    cout << "NNI assessed on: " << ((params.nni5) ? "5 branches" : "1 branch") << endl;
	}
	cout << "Phylogenetic likelihood library: " << (params.pll ? "Yes" : "No") << endl;
    cout << "Branch length optimization method: "
            << ((iqtree.optimize_by_newton) ? "Newton" : "Brent") << endl;
    cout << "Number of Newton-Raphson steps in NNI evaluation and branch length optimization: " << NNI_MAX_NR_STEP
            << " / " << PLL_NEWZPERCYCLE << endl;
    cout << "SSE instructions: "
            << ((iqtree.sse) ? "Yes" : "No") << endl;
	cout << endl;
}

void computeMLDist(Params& params, IQTree& iqtree, string &dist_file, double begin_time) {
	double longest_dist;
//	stringstream best_tree_string;
//	iqtree.printTree(best_tree_string, WT_BR_LEN + WT_TAXON_ID);
	cout << "Computing ML distances based on estimated model parameters...";
	double *ml_dist = NULL;
    double *ml_var = NULL;
    longest_dist = iqtree.computeDist(params, iqtree.aln, ml_dist, ml_var, dist_file);
	cout << " " << (getCPUTime() - begin_time) << " sec" << endl;
	if (longest_dist > MAX_GENETIC_DIST * 0.99) {
		outWarning("Some pairwise ML distances are too long (saturated)");
		//cout << "Some ML distances are too long, using old distances..." << endl;
	} //else
	{
		if ( !iqtree.dist_matrix ) {
	        iqtree.dist_matrix = new double[iqtree.aln->getNSeq() * iqtree.aln->getNSeq()];
		}
		if ( !iqtree.var_matrix ) {
	        iqtree.var_matrix = new double[iqtree.aln->getNSeq() * iqtree.aln->getNSeq()];
		}
		memmove(iqtree.dist_matrix, ml_dist,
                sizeof (double) * iqtree.aln->getNSeq() * iqtree.aln->getNSeq());
        memmove(iqtree.var_matrix, ml_var,
				sizeof(double) * iqtree.aln->getNSeq() * iqtree.aln->getNSeq());
	}
	delete[] ml_dist;
    delete[] ml_var;
}

void computeInitialDist(Params &params, IQTree &iqtree, string &dist_file) {
    double longest_dist;
	if (params.dist_file) {
		cout << "Reading distance matrix file " << params.dist_file << " ..." << endl;
	} else if (params.compute_jc_dist) {
		cout << "Computing Juke-Cantor distances..." << endl;
	} else if (params.compute_obs_dist) {
		cout << "Computing observed distances..." << endl;
	}

	if (params.compute_jc_dist || params.compute_obs_dist || params.partition_file) {
		longest_dist = iqtree.computeDist(params, iqtree.aln, iqtree.dist_matrix, iqtree.var_matrix, dist_file);
		checkZeroDist(iqtree.aln, iqtree.dist_matrix);
		if (longest_dist > MAX_GENETIC_DIST * 0.99) {
			outWarning("Some pairwise distances are too long (saturated)");
		}
    }

}

void initializeParams(Params &params, IQTree &iqtree, vector<ModelInfo> &model_info, ModelsBlock *models_block) {
//    iqtree.setCurScore(-DBL_MAX);
    bool test_only = params.model_name.find("ONLY") != string::npos;
    /* initialize substitution model */
    if (params.model_name.substr(0, 4) == "TEST") {
    	// TODO: check if necessary
//        if (iqtree.isSuperTree())
//            ((PhyloSuperTree*) &iqtree)->mapTrees();
        double start_cpu_time = getCPUTime();
        double start_real_time = getRealTime();
        ofstream fmodel;
        string fmodel_str = ((string)params.out_prefix + ".model"); 

        bool ok_model_file = false;
        if (!params.print_site_lh && !params.model_test_again) {
            ok_model_file = checkModelFile(fmodel_str, iqtree.isSuperTree(), model_info);
        }

        ok_model_file &= model_info.size() > 0;
        if (ok_model_file) {
            cout << "Reusing information from model file " << fmodel_str << endl;
            fmodel.open(fmodel_str.c_str(), ios::app);
            if (!fmodel.is_open())
                outError("cannot append to file ", fmodel_str);            
        } else {
            fmodel.open(fmodel_str.c_str());
            if (!fmodel.is_open())
                outError("cannot write to file ", fmodel_str);
            // print header
            SeqType seq_type = iqtree.aln->seq_type;
            if (iqtree.isSuperTree()) {
                fmodel << "Charset\t";
                seq_type = ((PhyloSuperTree*)&iqtree)->front()->aln->seq_type;
            }
            fmodel << "Model\tdf\tLnL\tTreeLen";
            if (seq_type == SEQ_BINARY)
                fmodel << "\t0\t1";
            else if (seq_type == SEQ_DNA)
                fmodel << "\tA-C\tA-G\tA-T\tC-G\tC-T\tG-T\tA\tC\tG\tT";
            fmodel << "\talpha\tpinv\tTree" << endl;
            model_info.clear();
        }
        fmodel.precision(4);
        fmodel << fixed;

        params.model_name = testModel(params, &iqtree, model_info, fmodel, models_block, "", true);
        fmodel.close();
        params.startCPUTime = start_cpu_time;
        params.start_real_time = start_real_time;
        cout << "CPU time for model selection: " << getCPUTime() - start_cpu_time << " seconds." << endl;
//        alignment = iqtree.aln;
        if (test_only) {
            params.min_iterations = 0;
        }
    }

    if (params.model_name == "WHTEST") {
        if (iqtree.aln->seq_type != SEQ_DNA)
            outError("Weiss & von Haeseler test of model homogeneity only works for DNA");
        params.model_name = "GTR+G";
    }

    assert(iqtree.aln);
    if (params.gbo_replicates)
        params.speed_conf = 1.0;

	// TODO: check if necessary
//    if (iqtree.isSuperTree())
//        ((PhyloSuperTree*) &iqtree)->mapTrees();

    // set parameter for the current tree
//    iqtree.setParams(params);
}


void pruneTaxa(Params &params, IQTree &iqtree, double *pattern_lh, NodeVector &pruned_taxa, StrVector &linked_name) {
	int num_low_support;
	double mytime;

	if (params.aLRT_threshold <= 100 && (params.aLRT_replicates > 0 || params.localbp_replicates > 0)) {
		mytime = getCPUTime();
		cout << "Testing tree branches by SH-like aLRT with " << params.aLRT_replicates << " replicates..." << endl;
		iqtree.setRootNode(params.root);
		double curScore =  iqtree.getCurScore();
		iqtree.computePatternLikelihood(pattern_lh, &curScore);
		num_low_support = iqtree.testAllBranches(params.aLRT_threshold, curScore,
				pattern_lh, params.aLRT_replicates, params.localbp_replicates, params.aLRT_test, params.aBayes_test);
		iqtree.printResultTree();
		cout << "  " << getCPUTime() - mytime << " sec." << endl;
		cout << num_low_support << " branches show low support values (<= " << params.aLRT_threshold << "%)" << endl;

		//tree.drawTree(cout);
		cout << "Collapsing stable clades..." << endl;
		iqtree.collapseStableClade(params.aLRT_threshold, pruned_taxa, linked_name, iqtree.dist_matrix);
		cout << pruned_taxa.size() << " taxa were pruned from stable clades" << endl;
	}

	if (!pruned_taxa.empty()) {
		cout << "Pruned alignment contains " << iqtree.aln->getNSeq()
				<< " sequences and " << iqtree.aln->getNSite() << " sites and "
				<< iqtree.aln->getNPattern() << " patterns" << endl;
		//tree.clearAllPartialLh();
		iqtree.initializeAllPartialLh();
		iqtree.clearAllPartialLH();
		iqtree.setCurScore(iqtree.optimizeAllBranches());
		//cout << "Log-likelihood	after reoptimizing model parameters: " << tree.curScore << endl;
		int nni_count, nni_steps;
		iqtree.setCurScore(iqtree.optimizeNNI(nni_count, nni_steps));
		cout << "Log-likelihood after optimizing partial tree: "
				<< iqtree.getCurScore() << endl;
	}

}

void restoreTaxa(IQTree &iqtree, double *saved_dist_mat, NodeVector &pruned_taxa, StrVector &linked_name) {
	if (!pruned_taxa.empty()) {
		cout << "Restoring full tree..." << endl;
		iqtree.restoreStableClade(iqtree.aln, pruned_taxa, linked_name);
		delete[] iqtree.dist_matrix;
		iqtree.dist_matrix = saved_dist_mat;
		iqtree.initializeAllPartialLh();
		iqtree.clearAllPartialLH();
		iqtree.setCurScore(iqtree.optimizeAllBranches());
		//cout << "Log-likelihood	after reoptimizing model parameters: " << tree.curScore << endl;
		int nni_count, nni_steps;
		iqtree.setCurScore(iqtree.optimizeNNI(nni_count, nni_steps));
		cout << "Log-likelihood	after reoptimizing full tree: " << iqtree.getCurScore() << endl;
		//iqtree.setBestScore(iqtree.getModelFactory()->optimizeParameters(params.fixed_branch_length, true, params.model_eps));

	}
}
void runApproximateBranchLengths(Params &params, IQTree &iqtree) {
    if (!params.fixed_branch_length && params.leastSquareBranch) {
        cout << endl << "Computing Least Square branch lengths..." << endl;
        iqtree.optimizeAllBranchesLS();
        iqtree.clearAllPartialLH();
        iqtree.setCurScore(iqtree.computeLikelihood());
        string filename = params.out_prefix;
        filename += ".lstree";
        iqtree.printTree(filename.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
        cout << "Logl of tree with LS branch lengths: " << iqtree.getCurScore() << endl;
        cout << "Tree with LS branch lengths written to " << filename << endl;
        if (params.print_branch_lengths) {
        	if (params.manuel_analytic_approx) {
        		cout << "Applying Manuel's analytic approximation.." << endl;
        		iqtree.approxAllBranches();
        	}
        	ofstream out;
        	filename = params.out_prefix;
        	filename += ".lsbrlen";
        	out.open(filename.c_str());
        	iqtree.printBranchLengths(out);
        	out.close();
        	cout << "LS Branch lengths written to " << filename << endl;
        }
        cout << "Total LS tree length: " << iqtree.treeLength() << endl;
    }

    if (params.pars_branch_length) {
    	cout << endl << "Computing parsimony branch lengths..." << endl;
    	iqtree.fixNegativeBranch(true);
    	iqtree.clearAllPartialLH();
        iqtree.setCurScore(iqtree.computeLikelihood());
        string filename = params.out_prefix;
        filename += ".mptree";
        iqtree.printTree(filename.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
        cout << "Logl of tree with MP branch lengths: " << iqtree.getCurScore() << endl;
        cout << "Tree with MP branch lengths written to " << filename << endl;
        if (params.print_branch_lengths) {
        	ofstream out;
        	filename = params.out_prefix;
        	filename += ".mpbrlen";
        	out.open(filename.c_str());
        	iqtree.printBranchLengths(out);
        	out.close();
        	cout << "MP Branch lengths written to " << filename << endl;
        }
        cout << "Total MP tree length: " << iqtree.treeLength() << endl;

    }

    if (params.bayes_branch_length) {
    	cout << endl << "Computing Bayesian branch lengths..." << endl;
    	iqtree.computeAllBayesianBranchLengths();
    	iqtree.clearAllPartialLH();
        iqtree.setCurScore(iqtree.computeLikelihood());
        string filename = params.out_prefix;
        filename += ".batree";
        iqtree.printTree(filename.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
        cout << "Logl of tree with Bayesian branch lengths: " << iqtree.getCurScore() << endl;
        cout << "Tree with Bayesian branch lengths written to " << filename << endl;
        if (params.print_branch_lengths) {
        	ofstream out;
        	filename = params.out_prefix;
        	filename += ".babrlen";
        	out.open(filename.c_str());
        	iqtree.printBranchLengths(out);
        	out.close();
        	cout << "Bayesian Branch lengths written to " << filename << endl;
        }
        cout << "Total Bayesian tree length: " << iqtree.treeLength() << endl;

    }

}

void printMiscInfo(Params &params, IQTree &iqtree, double *pattern_lh) {
	if (params.print_site_lh && !params.pll) {
		string site_lh_file = params.out_prefix;
		site_lh_file += ".sitelh";
		if (params.print_site_lh == WSL_SITE)
			printSiteLh(site_lh_file.c_str(), &iqtree, pattern_lh);
		else
			printSiteLhCategory(site_lh_file.c_str(), &iqtree, params.print_site_lh);
	}

    if (params.print_partition_lh && !iqtree.isSuperTree()) {
        outWarning("-wpl does not work with non-partition model");
        params.print_partition_lh = false;
    }
	if (params.print_partition_lh && !params.pll) {
        string part_lh_file = (string)params.out_prefix + ".partlh";
        printPartitionLh(part_lh_file.c_str(), &iqtree, pattern_lh);
	}

	if (params.print_site_prob && !params.pll) {
        printSiteProbCategory(((string)params.out_prefix + ".siteprob").c_str(), &iqtree, params.print_site_prob);
	}
    
    if (params.print_site_state_freq != WSF_NONE) {
		string site_freq_file = params.out_prefix;
		site_freq_file += ".sitesf";
        printSiteStateFreq(site_freq_file.c_str(), &iqtree);
    }

    if (params.print_trees_site_posterior) {
        cout << "Computing mixture posterior probabilities" << endl;
        IntVector pattern_cat;
        int num_mix = iqtree.computePatternCategories(&pattern_cat);
        cout << num_mix << " mixture components are necessary" << endl;
        string site_mix_file = (string)params.out_prefix + ".sitemix";
        ofstream out(site_mix_file.c_str());
        if (!out.is_open())
            outError("File " + site_mix_file + " could not be opened");
        out << "Ptn\tFreq\tNumMix" << endl;
        int ptn;
        for (ptn = 0; ptn < pattern_cat.size(); ptn++)
            out << ptn << "\t" << (int)iqtree.ptn_freq[ptn] << "\t" << pattern_cat[ptn] << endl;
        out.close();
        cout << "Pattern mixtures printed to " << site_mix_file << endl;
        
        site_mix_file = (string)params.out_prefix + ".sitemixall";
        out.open(site_mix_file.c_str());
        int ncat = iqtree.getRate()->getNRate();
        if (iqtree.getModel()->isMixture() && !iqtree.getModelFactory()->fused_mix_rate)
            ncat = iqtree.getModel()->getNMixtures();
        out << "Ptn\tFreq\tNumMix\tCat" << endl;
        
        int c;
        for (ptn = 0; ptn < iqtree.ptn_cat_mask.size(); ptn++) {
            int num_cat = popcount_lauradoux((unsigned*)&iqtree.ptn_cat_mask[ptn], 2);
            out << ptn << "\t" << (int)iqtree.ptn_freq[ptn] << "\t" << num_cat << "\t";
            for (c = 0; c < ncat; c++)
                if (iqtree.ptn_cat_mask[ptn] & ((uint64_t)1<<c))
                    out << "1";
                else
                    out << "0";
            out << endl;
        }
        out.close();
    }

	if (params.print_branch_lengths) {
    	if (params.manuel_analytic_approx) {
    		cout << "Applying Manuel's analytic approximation.." << endl;
    		iqtree.approxAllBranches();
    	}
		string brlen_file = params.out_prefix;
		brlen_file += ".brlen";
		ofstream out;
		out.open(brlen_file.c_str());
		iqtree.printBranchLengths(out);
		out.close();
		cout << "Branch lengths written to " << brlen_file << endl;
	}

	if (params.print_partition_info && iqtree.isSuperTree()) {
		string partition_info = params.out_prefix;
		partition_info += ".partinfo.nex";
		((PhyloSuperTree*)(&iqtree))->printPartition(partition_info.c_str());
		partition_info = (string)params.out_prefix + ".partitions";
		((PhyloSuperTree*)(&iqtree))->printPartitionRaxml(partition_info.c_str());
	}

	if (params.mvh_site_rate) {
		RateMeyerHaeseler *rate_mvh = new RateMeyerHaeseler(params.rate_file,
				&iqtree, params.rate_mh_type);
		cout << endl << "Computing site-specific rates by "
				<< rate_mvh->full_name << "..." << endl;
		rate_mvh->runIterativeProc(params, iqtree);
		cout << endl << "BEST SCORE FOUND : " << iqtree.candidateTrees.getBestScore()<< endl;
		string mhrate_file = params.out_prefix;
		mhrate_file += ".mhrate";
		iqtree.getRate()->writeSiteRates(mhrate_file.c_str());

		if (params.print_site_lh) {
			string site_lh_file = params.out_prefix;
			site_lh_file += ".mhsitelh";
			printSiteLh(site_lh_file.c_str(), &iqtree);
		}
	}

	if (params.print_site_rate) {
		string rate_file = params.out_prefix;
		rate_file += ".rate";
		iqtree.getRate()->writeSiteRates(rate_file.c_str());
		if (iqtree.isSuperTree()) {
			PhyloSuperTree *stree = (PhyloSuperTree*) &iqtree;
			int part = 0;
			try {
				ofstream out;
				out.exceptions(ios::failbit | ios::badbit);
				out.open(rate_file.c_str());
				for (PhyloSuperTree::iterator it = stree->begin(); it != stree->end(); it++, part++) {
					out << "SITE RATES FOR PARTITION " << stree->part_info[part].name << ":" << endl;
					(*it)->getRate()->writeSiteRates(out);
				}
				cout << "Site rates printed to " << rate_file << endl;
				out.close();
			} catch (ios::failure) {
				outError(ERR_WRITE_OUTPUT, rate_file);
			}
		}
	}

    if (params.fixed_branch_length == BRLEN_SCALE) {
        string filename = (string)params.out_prefix + ".blscale";
        iqtree.printTreeLengthScaling(filename.c_str());
        cout << "Scaled tree length and model parameters printed to " << filename << endl;
    }

}

void printFinalSearchInfo(Params &params, IQTree &iqtree, double search_cpu_time, double search_real_time) {
	cout << "Total tree length: " << iqtree.treeLength() << endl;

	if (iqtree.isSuperTree() && verbose_mode >= VB_MAX) {
		PhyloSuperTree *stree = (PhyloSuperTree*) &iqtree;
		cout << stree->evalNNIs << " NNIs evaluated from " << stree->totalNNIs << " all possible NNIs ( " <<
				(int)(((stree->evalNNIs+1.0)/(stree->totalNNIs+1.0))*100.0) << " %)" << endl;
		cout<<"Details for subtrees:"<<endl;
		for(int part = 0; part < stree->size(); part++){
			cout << part+1 <<". "<<stree->part_info[part].name<<": "<<stree->part_info[part].evalNNIs<<" ( "
				<< (int)(((stree->part_info[part].evalNNIs+1.0)/((stree->totalNNIs+1.0) / stree->size()))*100.0)
				<< " %)" << endl;
		}
	}

	params.run_time = (getCPUTime() - params.startCPUTime);
	cout << endl;
	cout << "Total number of iterations: " << iqtree.stop_rule.getCurIt() << endl;
//    cout << "Total number of partial likelihood vector computations: " << iqtree.num_partial_lh_computations << endl;
	cout << "CPU time used for tree search: " << search_cpu_time
			<< " sec (" << convert_time(search_cpu_time) << ")" << endl;
	cout << "Wall-clock time used for tree search: " << search_real_time
			<< " sec (" << convert_time(search_real_time) << ")" << endl;
	cout << "Total CPU time used: " << (double) params.run_time << " sec ("
			<< convert_time((double) params.run_time) << ")" << endl;
	cout << "Total wall-clock time used: "
			<< getRealTime() - params.start_real_time << " sec ("
			<< convert_time(getRealTime() - params.start_real_time) << ")" << endl;

}

void printSuboptimalTrees(IQTree& iqtree, Params& params, string suffix) {
	vector<string> trees = iqtree.candidateTrees.getTopTrees();
	ofstream treesOut((string(params.out_prefix) + suffix).c_str(),
			ofstream::out);
	for (vector<string>::iterator it = trees.begin(); it != trees.end(); it++) {
		treesOut << (*it);
		treesOut << endl;
	}
	treesOut.close();
}

/************************************************************
 *  MAIN TREE RECONSTRUCTION
 ***********************************************************/
void runTreeReconstruction(Params &params, string &original_model, IQTree &iqtree, vector<ModelInfo> &model_info) {

    string dist_file;
    params.startCPUTime = getCPUTime();
    params.start_real_time = getRealTime();

    // Make sure that no partial likelihood of IQ-TREE is initialized when PLL is used to save memory
    if (params.pll) {
        iqtree.deleteAllPartialLh();
    }

//    if (params.count_trees && pllTreeCounter == NULL)
//    	pllTreeCounter = new StringIntMap;

    // Temporary fix since PLL only supports DNA/Protein: switch to IQ-TREE parsimony kernel
    if (params.start_tree == STT_PLL_PARSIMONY) {
		if (iqtree.isSuperTree()) {
			PhyloSuperTree *stree = (PhyloSuperTree*)&iqtree;
			for (PhyloSuperTree::iterator it = stree->begin(); it != stree->end(); it++)
				if ((*it)->aln->seq_type != SEQ_DNA && (*it)->aln->seq_type != SEQ_PROTEIN)
					params.start_tree = STT_BIONJ;
		} else if (iqtree.aln->seq_type != SEQ_DNA && iqtree.aln->seq_type != SEQ_PROTEIN)
			params.start_tree = STT_PARSIMONY;
    }

    /***************** Initialization for PLL and sNNI ******************/
    if (params.start_tree == STT_PLL_PARSIMONY || params.start_tree == STT_RANDOM_TREE || params.pll) {
        /* Initialized all data structure for PLL*/
    	iqtree.initializePLL(params);
    }


    /********************* Compute pairwise distances *******************/
    if (params.start_tree == STT_BIONJ || params.iqp || params.leastSquareBranch) {
    	computeInitialDist(params, iqtree, dist_file);
    }

    /******************** Pass the parameter object params to IQTree *******************/
    iqtree.setParams(&params);

    /********************** Create an initial tree **********************/
    iqtree.computeInitialTree(dist_file, params.SSE);
    
    //*** FOR TUNG: This is wrong! a NULL root was already treated correctly
//    if (params.root == NULL) {
//    	params.root = iqtree.aln->getSeqName(0).c_str();
//    	iqtree.setRootNode(params.root);
//    }
   	iqtree.setRootNode(params.root);

    /*************** SET UP PARAMETERS and model testing ****************/

   	// FOR TUNG: swapping the order cause bug for -m TESTLINK
//    iqtree.initSettings(params);

	ModelsBlock *models_block = readModelsDefinition(params);

    initializeParams(params, iqtree, model_info, models_block);

    iqtree.restoreCheckpoint();
    iqtree.initSettings(params);

    /*********************** INITIAL MODEL OPTIMIZATION *****************/

    iqtree.initializeModel(params, models_block);

    delete models_block;

    // UpperBounds analysis. Here, to analyse the initial tree without any tree search or optimization
    if (params.upper_bound) {
    	iqtree.setCurScore(iqtree.computeLikelihood());
    	cout<<iqtree.getCurScore()<<endl;
    	UpperBounds(&params, iqtree.aln, &iqtree);
    	exit(0);
	}

    // degree of freedom
    cout << endl;
    if (verbose_mode >= VB_MED) {
    	cout << "ML-TREE SEARCH START WITH THE FOLLOWING PARAMETERS:" << endl;
        int model_df = iqtree.getModelFactory()->getNParameters();
    	printAnalysisInfo(model_df, iqtree, params);
    }

    if (!params.pll) {
        uint64_t mem_size = iqtree.getMemoryRequired();
        uint64_t total_mem = getMemorySize();
        if (mem_size >= total_mem) {
            if (params.lh_mem_save == LM_DETECT) {
                // switch to memory saving technique that reduces memory requirement to 1/3
                params.lh_mem_save = LM_PER_NODE;
                mem_size = iqtree.getMemoryRequired();
            }
        }
//#if defined __APPLE__ || defined __MACH__
        cout << "NOTE: " << (mem_size / 1024) / 1024 << " MB RAM is required!" << endl;
//#else
//        cout << "NOTE: " << ((double) mem_size / 1000.0) / 1000 << " MB RAM is required!" << endl;
//#endif
        if (mem_size >= total_mem) {
            outError("Memory required exceeds your computer RAM size!");
        }
#ifdef BINARY32
        if (mem_size >= 2000000000) {
            outError("Memory required exceeds 2GB limit of 32-bit executable");
        }
#endif
        int max_procs = countPhysicalCPUCores();
        if (mem_size * max_procs > total_mem * params.num_threads) {
            outWarning("Memory required per CPU-core (" + convertDoubleToString((double)mem_size/params.num_threads/1024/1024/1024)+
            " GB) is higher than your computer RAM per CPU-core ("+convertIntToString(total_mem/max_procs/1024/1024/1024)+
            " GB), thus multiple runs may exceed RAM!");
        }
    }

    iqtree.initializeAllPartialLh();
	double initEpsilon = params.min_iterations == 0 ? params.modeps : (params.modeps*10);

	string initTree;

	if (iqtree.getRate()->name.find("+I+G") != string::npos) {
		if (params.alpha_invar_file != NULL) { // COMPUTE TREE LIKELIHOOD BASED ON THE INPUT ALPHA AND P_INVAR VALUE
			computeLoglFromUserInputGAMMAInvar(params, iqtree);
			exit(0);
		}

		if (params.exh_ai) {
			exhaustiveSearchGAMMAInvar(params, iqtree);
			exit(0);
		}

	}

    // Optimize model parameters and branch lengths using ML for the initial tree
	iqtree.clearAllPartialLH();
    iqtree.getModelFactory()->restoreCheckpoint();
    if (iqtree.getCheckpoint()->getBool("finishedModelInit")) {
        // model optimization already done: ignore this step
        if (!iqtree.candidateTrees.empty())
            iqtree.readTreeString(iqtree.candidateTrees.getTopTrees(1)[0]);
        iqtree.setCurScore(iqtree.computeLikelihood());
        initTree = iqtree.getTreeString();
        cout << "CHECKPOINT: Model parameters restored, LogL: " << iqtree.getCurScore() << endl;
    } else {
        initTree = iqtree.optimizeModelParameters(true, initEpsilon);
        iqtree.saveCheckpoint();
        iqtree.getModelFactory()->saveCheckpoint();
        iqtree.getCheckpoint()->putBool("finishedModelInit", true);
        iqtree.getCheckpoint()->dump();
    }

    if (params.lmap_num_quartets >= 0) {
        cout << endl << "Performing likelihood mapping with ";
        if (params.lmap_num_quartets > 0)
            cout << params.lmap_num_quartets;
        else
            cout << "all";
        cout << " quartets..." << endl;
        double lkmap_time = getRealTime();
        iqtree.doLikelihoodMapping();
        cout << "Likelihood mapping needed " << getRealTime()-lkmap_time << " seconds" << endl << endl;
    }
    
    bool finishedCandidateSet = iqtree.getCheckpoint()->getBool("finishedCandidateSet");
    bool finishedInitTree = iqtree.getCheckpoint()->getBool("finishedInitTree");

    // now overwrite with random tree
    if (params.start_tree == STT_RANDOM_TREE && !finishedInitTree) {
        cout << "Generate random initial Yule-Harding tree..." << endl;
        iqtree.generateRandomTree(YULE_HARDING);
        iqtree.wrapperFixNegativeBranch(true);
        iqtree.initializeAllPartialLh();
        initTree = iqtree.optimizeBranches(2);
        cout << "Log-likelihood of random tree: " << iqtree.getCurScore() << endl;
    }

    /****************** NOW PERFORM MAXIMUM LIKELIHOOD TREE RECONSTRUCTION ******************/

    // Update best tree
    if (!finishedInitTree)
        iqtree.candidateTrees.update(initTree, iqtree.getCurScore());

    if (params.min_iterations > 0) {
        if (!iqtree.isBifurcating())
            outError("Tree search does not work with initial multifurcating tree. Please specify `-n 0` to avoid this.");
        cout << "--------------------------------------------------------------------" << endl;
        cout << "|             INITIALIZING CANDIDATE TREE SET                      |" << endl;
        cout << "--------------------------------------------------------------------" << endl;
    }

    // Compute maximum likelihood distance
    // ML distance is only needed for IQP
//    if ( params.start_tree != STT_BIONJ && ((params.snni && !params.iqp) || params.min_iterations == 0)) {
//        params.compute_ml_dist = false;
//    }
    if (params.min_iterations <= 1 && params.start_tree != STT_BIONJ)
        params.compute_ml_dist = false;
    
    if ((params.user_file || params.start_tree == STT_RANDOM_TREE) && params.snni && !params.iqp) {
        params.compute_ml_dist = false;
    }
    
    if (params.constraint_tree_file)
        params.compute_ml_dist = false;

//    if ( params.user_file && params.min_iterations == 0) {
//        params.compute_ml_dist = false;
//    }

    if (!finishedInitTree && ((!params.dist_file && params.compute_ml_dist) || params.leastSquareBranch)) {
        computeMLDist(params, iqtree, dist_file, getCPUTime());
        if (!params.user_file && params.start_tree != STT_RANDOM_TREE) {
            // NEW 2015-08-10: always compute BIONJ tree into the candidate set
            iqtree.resetCurScore();
            double start_bionj = getRealTime();
            iqtree.computeBioNJ(params, iqtree.aln, dist_file);
            cout << getRealTime() - start_bionj << " seconds" << endl;
            if (iqtree.isSuperTree())
                iqtree.wrapperFixNegativeBranch(true);
            else
                iqtree.wrapperFixNegativeBranch(false);
            if (params.start_tree == STT_BIONJ) {
                initTree = iqtree.optimizeModelParameters(params.min_iterations==0, initEpsilon);
            } else {
                initTree = iqtree.optimizeBranches();
            }
            cout << "Log-likelihood of BIONJ tree: " << iqtree.getCurScore() << endl;
            iqtree.candidateTrees.update(initTree, iqtree.getCurScore());
        }
    }

//    iqtree.saveCheckpoint();

	double cputime_search_start = getCPUTime();
    double realtime_search_start = getRealTime();

    if (params.min_iterations > 0 && !finishedCandidateSet) {
        double initTime = getCPUTime();

//        if (!params.user_file && (params.start_tree == STT_PARSIMONY || params.start_tree == STT_PLL_PARSIMONY)) 
//        {
        	iqtree.initCandidateTreeSet(params.numInitTrees - iqtree.candidateTrees.size(), params.numNNITrees);
        	assert(iqtree.candidateTrees.size() != 0);
        	cout << "Finish initializing candidate tree set. ";
        	cout << "Number of distinct locally optimal trees: " << iqtree.candidateTrees.size() << endl;
        	if (params.write_local_optimal_trees) {
        		printSuboptimalTrees(iqtree, params, ".init_suboptimal_trees");
        	}
//        }
        cout << "Current best tree score: " << iqtree.candidateTrees.getBestScore() << " / CPU time: "
                << getCPUTime() - initTime << endl;
	}

    if (finishedCandidateSet) {
        cout << "CHECKPOINT: Candidate tree set restored, best LogL: " << iqtree.candidateTrees.getBestScore() << endl;
    } else {
        iqtree.saveCheckpoint();
        iqtree.getCheckpoint()->putBool("finishedCandidateSet", true);
        iqtree.getCheckpoint()->dump(true);
    }

    if (params.leastSquareNNI) {
    	iqtree.computeSubtreeDists();
    }
    /* TUNG: what happens if params.root is not set? This is usually the case.
     * I added code to ininialize the root above.
     */
    //iqtree.setRootNode(params.root); // Important for NNI below

	if (original_model == "WHTEST") {
		cout << endl << "Testing model homogeneity by Weiss & von Haeseler (2003)..." << endl;
		WHTest(params, iqtree);
	}

	NodeVector pruned_taxa;
	StrVector linked_name;
	double *saved_dist_mat = iqtree.dist_matrix;
	double *pattern_lh;

	pattern_lh = new double[iqtree.getAlnNPattern()];

	// prune stable taxa
	pruneTaxa(params, iqtree, pattern_lh, pruned_taxa, linked_name);

	if (params.min_iterations > 1) {
		iqtree.readTreeString(iqtree.candidateTrees.getTopTrees()[0]);
		iqtree.doTreeSearch();
		iqtree.setAlignment(iqtree.aln);
        cout << "TREE SEARCH COMPLETED AFTER " << iqtree.stop_rule.getCurIt() << " ITERATIONS" 
            << " / Time: " << convert_time(getRealTime() - params.start_real_time) << endl << endl;
	} else {
		/* do SPR with likelihood function */
		if (params.tree_spr) {
			//tree.optimizeSPRBranches();
			cout << "Doing SPR Search" << endl;
			cout << "Start tree.optimizeSPR()" << endl;
			double spr_score = iqtree.optimizeSPR();
			cout << "Finish tree.optimizeSPR()" << endl;
			//double spr_score = tree.optimizeSPR(tree.curScore, (PhyloNode*) tree.root->neighbors[0]->node);
			if (spr_score <= iqtree.getCurScore()) {
				cout << "SPR search did not found any better tree" << endl;
			}
		}
	}

	// restore pruned taxa
	restoreTaxa(iqtree, saved_dist_mat, pruned_taxa, linked_name);

	double search_cpu_time = getCPUTime() - cputime_search_start;
	double search_real_time = getRealTime() - realtime_search_start;

    // COMMENT THIS OUT BECAUSE IT DELETES ALL BRANCH LENGTHS OF SUBTREES!
//	if (iqtree.isSuperTree())
//			((PhyloSuperTree*) &iqtree)->mapTrees();

	if (params.snni && params.min_iterations && verbose_mode >= VB_MED) {
		cout << "Log-likelihoods of best " << params.popSize << " trees: " << endl;
		iqtree.printBestScores(params.popSize);
		cout << endl;
	}

	if (params.min_iterations) {
		iqtree.readTreeString(iqtree.candidateTrees.getBestTrees()[0]);
        iqtree.initializeAllPartialLh();
        iqtree.clearAllPartialLH();
        cout << "--------------------------------------------------------------------" << endl;
        cout << "|                    FINALIZING TREE SEARCH                        |" << endl;
        cout << "--------------------------------------------------------------------" << endl;
        
        if (iqtree.getCheckpoint()->getBool("finishedModelFinal")) {
            iqtree.setCurScore(iqtree.computeLikelihood());
            cout << "CHECKPOINT: Final model parameters restored" << endl;
        } else {
            cout << "Performs final model parameters optimization" << endl;
            string tree;
            tree = iqtree.optimizeModelParameters(true);
            iqtree.candidateTrees.update(tree, iqtree.getCurScore(), true);
            iqtree.getCheckpoint()->putBool("finishedModelFinal", true);
            iqtree.saveCheckpoint();
        }
        
    }

	if (iqtree.isSuperTree())
		((PhyloSuperTree*) &iqtree)->computeBranchLengths();

	cout << "BEST SCORE FOUND : " << iqtree.getCurScore() << endl;

	if (params.write_local_optimal_trees) {
		printSuboptimalTrees(iqtree, params, ".suboptimal_trees");
	}

	if (params.pll)
		iqtree.inputModelPLL2IQTree();

	/* root the tree at the first sequence */
	iqtree.root = iqtree.findLeafName(iqtree.aln->getSeqName(0));
	assert(iqtree.root);


	if (!params.pll) {
	    iqtree.computeLikelihood(pattern_lh);
	    // compute logl variance
	    iqtree.logl_variance = iqtree.computeLogLVariance();
	}

	printMiscInfo(params, iqtree, pattern_lh);

	/****** perform SH-aLRT test ******************/
	if ((params.aLRT_replicates > 0 || params.localbp_replicates > 0 || params.aLRT_test || params.aBayes_test) && !params.pll) {
		double mytime = getCPUTime();
		params.aLRT_replicates = max(params.aLRT_replicates, params.localbp_replicates);
        cout << endl;
        if (params.aLRT_replicates > 0) 
            cout << "Testing tree branches by SH-like aLRT with "
				<< params.aLRT_replicates << " replicates..." << endl;
        if (params.localbp_replicates)
            cout << "Testing tree branches by local-BP test with " << params.localbp_replicates << " replicates..." << endl;
        if (params.aLRT_test)
            cout << "Testing tree branches by aLRT parametric test..." << endl;
        if (params.aBayes_test)
            cout << "Testing tree branches by aBayes parametric test..." << endl;
		iqtree.setRootNode(params.root);
        if (iqtree.isBifurcating()) {
            iqtree.testAllBranches(params.aLRT_threshold, iqtree.getCurScore(),
                    pattern_lh, params.aLRT_replicates, params.localbp_replicates, params.aLRT_test, params.aBayes_test);
            cout << "CPU Time used:  " << getCPUTime() - mytime << " sec." << endl;
        } else {
            outWarning("Tree is multifurcating and such test is not applicable");
            params.aLRT_replicates = params.localbp_replicates = params.aLRT_test = params.aBayes_test = 0;
        }
	}

	if (params.gbo_replicates > 0) {
		if (!params.online_bootstrap)
			outError("Obsolete feature");
//			runGuidedBootstrap(params, iqtree.aln, iqtree);
		else
			iqtree.summarizeBootstrap(params);
	}

	printFinalSearchInfo(params, iqtree, search_cpu_time, search_real_time);

	// BUG FIX: readTreeString(bestTreeString) not needed before this line
	iqtree.printResultTree();

	if(params.upper_bound_NNI){
		string out_file_UB = params.out_prefix;
		out_file_UB += ".UB.NNI.main";
		ofstream out_UB;
		out_UB.exceptions(ios::failbit | ios::badbit);
		out_UB.open((char*)out_file_UB.c_str(),std::ofstream::out | std::ofstream::app);
		out_UB<<iqtree.leafNum<<"\t"<<iqtree.aln->getNSite()<<"\t"<<iqtree.params->upper_bound_frac<<"\t"
				  <<iqtree.skippedNNIub<<"\t"<< iqtree.totalNNIub<<"\t"<<iqtree.candidateTrees.getBestScore() <<endl;
					//iqtree.minUB << "\t" << iqtree.meanUB/iqtree.skippedNNIub << "\t" << iqtree.maxUB << endl;
		out_UB.close();
		}

	if (params.out_file)
		iqtree.printTree(params.out_file);

	delete[] pattern_lh;

	runApproximateBranchLengths(params, iqtree);

}

void computeLoglFromUserInputGAMMAInvar(Params &params, IQTree &iqtree) {
	RateHeterogeneity *site_rates = iqtree.getRate();
	site_rates->setFixPInvar(true);
	site_rates->setFixGammaShape(true);
	vector<double> alphas, p_invars, logl;
	ifstream aiFile;
	aiFile.open(params.alpha_invar_file, ios_base::in);
	if (aiFile.good()) {
		double alpha, p_invar;
		while (aiFile >> alpha >> p_invar) {
			alphas.push_back(alpha);
			p_invars.push_back(p_invar);
		}
		aiFile.close();
		cout << "Computing tree logl based on the alpha and p_invar values in " << params.alpha_invar_file << " ..." <<
		endl;
	} else {
		stringstream errMsg;
		errMsg << "Could not find file: " << params.alpha_invar_file;
		outError(errMsg.str().c_str());
	}
	string aiResultsFileName = string(params.out_prefix) + "_" + string(params.alpha_invar_file) + ".results";
	ofstream aiFileResults;
	aiFileResults.open(aiResultsFileName.c_str());
	aiFileResults << fixed;
	aiFileResults.precision(4);
	DoubleVector lenvec;
	aiFileResults << "Alpha P_Invar Logl TreeLength\n";
	for (int i = 0; i < alphas.size(); i++) {
		iqtree.saveBranchLengths(lenvec);
		aiFileResults << alphas.at(i) << " " << p_invars.at(i) << " ";
		site_rates->setGammaShape(alphas.at(i));
		site_rates->setPInvar(p_invars.at(i));
		iqtree.clearAllPartialLH();
		double lh = iqtree.getModelFactory()->optimizeParameters(params.fixed_branch_length, false, 0.001);
		aiFileResults << lh << " " << iqtree.treeLength() << "\n";
		iqtree.restoreBranchLengths(lenvec);
	}
	aiFileResults.close();
	cout << "Results were written to: " << aiResultsFileName << endl;
	cout << "Wall clock time used: " << getRealTime() - params.start_real_time << endl;
}

void searchGAMMAInvarByRestarting(IQTree &iqtree) {
    if (!Params::getInstance().fixed_branch_length)
		iqtree.setCurScore(iqtree.optimizeAllBranches(1));
	else
		iqtree.setCurScore(iqtree.computeLikelihood());
	RateHeterogeneity* site_rates = (iqtree.getRate());
	double values[] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
	vector<double> initAlphas;
	if (Params::getInstance().randomAlpha) {
		while (initAlphas.size() < 10) {
			double initAlpha = random_double();
			initAlphas.push_back(initAlpha + MIN_GAMMA_SHAPE*2);
		}
	} else {
		initAlphas.assign(values, values+10);
	}
	double bestLogl = iqtree.getCurScore();
	double bestAlpha = 0.0;
	double bestPInvar = 0.0;
	double initPInvar = iqtree.getRate()->getPInvar();

    /* Back up branch lengths and substitutional rates */
	DoubleVector lenvec;
	DoubleVector bestLens;
	iqtree.saveBranchLengths(lenvec);
    int numRateEntries = iqtree.getModel()->getNumRateEntries();
    double *rates = new double[numRateEntries];
    double *bestRates = new double[numRateEntries];
    iqtree.getModel()->getRateMatrix(rates);
    int numStates = iqtree.aln->num_states;
    double *state_freqs = new double[numStates];
    iqtree.getModel()->getStateFrequency(state_freqs);
    double *bestStateFreqs =  new double[numStates];

    for (int i = 0; i < 10; i++) {
		cout << endl;
		cout << "Testing alpha: " << initAlphas[i] << endl;
        // Initialize model parameters
        iqtree.restoreBranchLengths(lenvec);
        ((ModelGTR*) iqtree.getModel())->setRateMatrix(rates);
        ((ModelGTR*) iqtree.getModel())->setStateFrequency(state_freqs);
        iqtree.getModel()->decomposeRateMatrix();
        site_rates->setGammaShape(initAlphas[i]);
		site_rates->setPInvar(initPInvar);
		iqtree.clearAllPartialLH();
		iqtree.optimizeModelParameters(verbose_mode >= VB_MED, Params::getInstance().testAlphaEps);
        double estAlpha = iqtree.getRate()->getGammaShape();
        double estPInv = iqtree.getRate()->getPInvar();
        double logl = iqtree.getCurScore();
		cout << "Est. alpha: " << estAlpha << " / Est. pinv: " << estPInv
        << " / Logl: " << logl << endl;

		if (iqtree.getCurScore() > bestLogl) {
			bestLogl = logl;
			bestAlpha = estAlpha;
			bestPInvar = estPInv;
			bestLens.clear();
			iqtree.saveBranchLengths(bestLens);
            iqtree.getModel()->getRateMatrix(bestRates);
            iqtree.getModel()->getStateFrequency(bestStateFreqs);
        }
    }
	site_rates->setGammaShape(bestAlpha);
	site_rates->setFixGammaShape(false);
	site_rates->setPInvar(bestPInvar);
	site_rates->setFixPInvar(false);
    ((ModelGTR*) iqtree.getModel())->setRateMatrix(bestRates);
    ((ModelGTR*) iqtree.getModel())->setStateFrequency(bestStateFreqs);
	iqtree.restoreBranchLengths(bestLens);
    iqtree.getModel()->decomposeRateMatrix();
	iqtree.clearAllPartialLH();
    iqtree.setCurScore(iqtree.computeLikelihood());
    cout << endl;
    cout << "Best initial alpha: " << bestAlpha << " / initial pinv: " << bestPInvar << " / ";
    cout << "Logl: " << iqtree.getCurScore() << endl;

    delete [] rates;
    delete [] state_freqs;
    delete [] bestRates;
    delete [] bestStateFreqs;
}

// Test alpha fom 0.1 to 15 and p_invar from 0.1 to 0.99, stepsize = 0.01
void exhaustiveSearchGAMMAInvar(Params &params, IQTree &iqtree) {
	double alphaMin = 0.01;
	double alphaMax = 2.00;
	double p_invarMin = 0.01;
	double p_invarMax = 1.00;
	double stepSize = 0.01;
	int numAlpha = (int) floor((alphaMax - alphaMin)/stepSize);
	int numInvar = (int) floor((p_invarMax - p_invarMin)/stepSize);

	cout << "EVALUATING COMBINATIONS OF " << numAlpha << " ALPHAS AND " << numInvar << " P_INVARS ... " << endl;

	vector<string> results;
	results.reserve((unsigned long) (numAlpha * numInvar));
	DoubleVector lenvec;
	iqtree.saveBranchLengths(lenvec);

	RateHeterogeneity* site_rates = (iqtree.getRate());
	site_rates->setFixPInvar(true);
	site_rates->setFixGammaShape(true);

    for (double alpha = alphaMin; alpha < alphaMax; alpha = alpha + stepSize) {
        for (double p_invar = p_invarMin; p_invar < p_invarMax; p_invar = p_invar + stepSize) {
            site_rates->setGammaShape(alpha);
            site_rates->setPInvar(p_invar);
            iqtree.clearAllPartialLH();
            double lh = iqtree.getModelFactory()->optimizeParameters(params.fixed_branch_length, false, 0.001);
            stringstream ss;
            ss << fixed << setprecision(2) << alpha << " " << p_invar << " " << lh << " " << iqtree.treeLength();
            //cout << ss.str() << endl;
            results.push_back(ss.str());
            iqtree.restoreBranchLengths(lenvec);
        }
    }
	string aiResultsFileName = string(params.out_prefix) + ".ai_results";
	ofstream aiFileResults;
	aiFileResults.open(aiResultsFileName.c_str());
	aiFileResults << fixed;
	aiFileResults.precision(4);
	aiFileResults << "alpha p_invar logl tree_len\n";
	for (vector<string>::iterator it = results.begin(); it != results.end(); it++) {
				aiFileResults << (*it) << endl;
			}
	aiFileResults.close();
	cout << "Results were written to: " << aiResultsFileName << endl;
	cout << "Wall clock time used: " << getRealTime() - params.start_real_time << endl;
}


/**********************************************************
 * STANDARD NON-PARAMETRIC BOOTSTRAP
 ***********************************************************/
void runStandardBootstrap(Params &params, string &original_model, Alignment *alignment, IQTree *tree) {
	vector<ModelInfo> *model_info = new vector<ModelInfo>;
	StrVector removed_seqs, twin_seqs;

	// turn off all branch tests
	int saved_aLRT_replicates = params.aLRT_replicates;
    int saved_localbp_replicates = params.localbp_replicates;
    bool saved_aLRT_test = params.aLRT_test;
    bool saved_aBayes_test = params.aBayes_test;
	params.aLRT_replicates = 0;
    params.localbp_replicates = 0;
    params.aLRT_test = false;
    params.aBayes_test = false;
    
    if (params.suppress_output_flags & OUT_TREEFILE)
        outError("Suppress .treefile not allowed for standard bootstrap");
	string treefile_name = params.out_prefix;
	treefile_name += ".treefile";
	string boottrees_name = params.out_prefix;
	boottrees_name += ".boottrees";
	string bootaln_name = params.out_prefix;
	bootaln_name += ".bootaln";
	string bootlh_name = params.out_prefix;
	bootlh_name += ".bootlh";
    int bootSample = 0;
    if (tree->getCheckpoint()->get("bootSample", bootSample)) {
        cout << "CHECKPOINT: " << bootSample << " bootstrap analyses restored" << endl;
    } else {
        // first empty the boottrees file
        try {
            ofstream tree_out;
            tree_out.exceptions(ios::failbit | ios::badbit);
            tree_out.open(boottrees_name.c_str());
            tree_out.close();
        } catch (ios::failure) {
            outError(ERR_WRITE_OUTPUT, boottrees_name);
        }

        // empty the bootaln file
        if (params.print_bootaln)
        try {
            ofstream tree_out;
            tree_out.exceptions(ios::failbit | ios::badbit);
            tree_out.open(bootaln_name.c_str());
            tree_out.close();
        } catch (ios::failure) {
            outError(ERR_WRITE_OUTPUT, bootaln_name);
        }
    }
    
	double start_time = getCPUTime();

    
    
	// do bootstrap analysis
	for (int sample = bootSample; sample < params.num_bootstrap_samples; sample++) {
		cout << endl << "===> START BOOTSTRAP REPLICATE NUMBER "
				<< sample + 1 << endl << endl;

        // 2015-12-17: initialize random stream for creating bootstrap samples
        // mainly so that checkpointing does not need to save bootstrap samples
        int *saved_randstream = randstream;
        init_random(params.ran_seed + sample);

		Alignment* bootstrap_alignment;
		cout << "Creating bootstrap alignment (seed: " << params.ran_seed+sample << ")..." << endl;
		if (alignment->isSuperAlignment())
			bootstrap_alignment = new SuperAlignment;
		else
			bootstrap_alignment = new Alignment;
		bootstrap_alignment->createBootstrapAlignment(alignment, NULL, params.bootstrap_spec);

        // restore randstream
        finish_random();
        randstream = saved_randstream;

		if (params.print_tree_lh) {
			double prob;
			bootstrap_alignment->multinomialProb(*alignment, prob);
			ofstream boot_lh;
			if (sample == 0)
				boot_lh.open(bootlh_name.c_str());
			else
				boot_lh.open(bootlh_name.c_str(), ios_base::out | ios_base::app);
			boot_lh << "0\t" << prob << endl;
			boot_lh.close();
		}
		IQTree *boot_tree;
		if (alignment->isSuperAlignment()){
			if(params.partition_type){
				boot_tree = new PhyloSuperTreePlen((SuperAlignment*) bootstrap_alignment, (PhyloSuperTree*) tree);
			} else {
				boot_tree = new PhyloSuperTree((SuperAlignment*) bootstrap_alignment, (PhyloSuperTree*) tree);
			}
		} else
			boot_tree = new IQTree(bootstrap_alignment);
		if (params.print_bootaln)
			bootstrap_alignment->printPhylip(bootaln_name.c_str(), true);

        // set checkpoint
        boot_tree->setCheckpoint(tree->getCheckpoint());
        boot_tree->num_precision = tree->num_precision;

		runTreeReconstruction(params, original_model, *boot_tree, *model_info);
		// read in the output tree file
        stringstream ss;
        boot_tree->printTree(ss);
//		try {
//			ifstream tree_in;
//			tree_in.exceptions(ios::failbit | ios::badbit);
//			tree_in.open(treefile_name.c_str());
//			tree_in >> tree_str;
//			tree_in.close();
//		} catch (ios::failure) {
//			outError(ERR_READ_INPUT, treefile_name);
//		}
		// write the tree into .boottrees file
		try {
			ofstream tree_out;
			tree_out.exceptions(ios::failbit | ios::badbit);
			tree_out.open(boottrees_name.c_str(), ios_base::out | ios_base::app);
			tree_out << ss.str() << endl;
			tree_out.close();
		} catch (ios::failure) {
			outError(ERR_WRITE_OUTPUT, boottrees_name);
		}
		// fix bug: set the model for original tree after testing
		if (original_model.substr(0,4) == "TEST" && tree->isSuperTree()) {
			PhyloSuperTree *stree = ((PhyloSuperTree*)tree);
			stree->part_info =  ((PhyloSuperTree*)boot_tree)->part_info;
//			for (int i = 0; i < ((PhyloSuperTree*)tree)->part_info.size(); i++)
//				((PhyloSuperTree*)tree)->part_info[i].model_name = ((PhyloSuperTree*)boot_tree)->part_info[i].model_name;
		}
		if (params.num_bootstrap_samples == 1)
			reportPhyloAnalysis(params, original_model, *boot_tree, *model_info);
		// WHY was the following line missing, which caused memory leak?
		bootstrap_alignment = boot_tree->aln;
		delete boot_tree;
		// fix bug: bootstrap_alignment might be changed
		delete bootstrap_alignment;
        
        // clear all checkpointed information
        Checkpoint *newCheckpoint = new Checkpoint;
        tree->getCheckpoint()->getSubCheckpoint(newCheckpoint, "iqtree");
        tree->getCheckpoint()->clear();
        tree->getCheckpoint()->insert(newCheckpoint->begin(), newCheckpoint->end());
        tree->getCheckpoint()->put("bootSample", sample+1);
        tree->getCheckpoint()->putBool("finished", false);
        tree->getCheckpoint()->dump(true);
        delete newCheckpoint;
        
	}


	if (params.consensus_type == CT_CONSENSUS_TREE) {

		cout << endl << "===> COMPUTE CONSENSUS TREE FROM "
				<< params.num_bootstrap_samples << " BOOTSTRAP TREES" << endl << endl;
		computeConsensusTree(boottrees_name.c_str(), 0, 1e6, -1,
				params.split_threshold, NULL, params.out_prefix, NULL, &params);
	}

	if (params.compute_ml_tree) {
		cout << endl << "===> START ANALYSIS ON THE ORIGINAL ALIGNMENT" << endl << endl;
        // restore branch tests
		params.aLRT_replicates = saved_aLRT_replicates;
        params.localbp_replicates = saved_localbp_replicates;
        params.aLRT_test = saved_aLRT_test;
        params.aBayes_test = saved_aBayes_test;
        
		runTreeReconstruction(params, original_model, *tree, *model_info);

		cout << endl << "===> ASSIGN BOOTSTRAP SUPPORTS TO THE TREE FROM ORIGINAL ALIGNMENT" << endl << endl;
		MExtTree ext_tree;
		assignBootstrapSupport(boottrees_name.c_str(), 0, 1e6,
				treefile_name.c_str(), false, treefile_name.c_str(),
				params.out_prefix, ext_tree, NULL, &params);
		tree->copyTree(&ext_tree);
		reportPhyloAnalysis(params, original_model, *tree, *model_info);
	} else if (params.consensus_type == CT_CONSENSUS_TREE) {
		int mi = params.min_iterations;
		STOP_CONDITION sc = params.stop_condition;
		params.min_iterations = 0;
		params.stop_condition = SC_FIXED_ITERATION;
		runTreeReconstruction(params, original_model, *tree, *model_info);
		params.min_iterations = mi;
		params.stop_condition = sc;
		tree->stop_rule.initialize(params);
		reportPhyloAnalysis(params, original_model, *tree, *model_info);
	} else
		cout << endl;

	cout << "Total CPU time for bootstrap: " << (getCPUTime() - start_time) << " seconds." << endl << endl;
	cout << "Non-parametric bootstrap results written to:" << endl;
	if (params.print_bootaln)
		cout << "  Bootstrap alignments:     " << params.out_prefix << ".bootaln" << endl;
	cout << "  Bootstrap trees:          " << params.out_prefix << ".boottrees" << endl;
	if (params.consensus_type == CT_CONSENSUS_TREE)
		cout << "  Consensus tree:           " << params.out_prefix << ".contree" << endl;
	cout << endl;
    
    delete model_info;
}

void convertAlignment(Params &params, IQTree *iqtree) {
	Alignment *alignment = iqtree->aln;
	if (params.num_bootstrap_samples || params.print_bootaln) {
		// create bootstrap alignment
		Alignment* bootstrap_alignment;
		cout << "Creating bootstrap alignment..." << endl;
		if (alignment->isSuperAlignment())
			bootstrap_alignment = new SuperAlignment;
		else
			bootstrap_alignment = new Alignment;
		bootstrap_alignment->createBootstrapAlignment(alignment, NULL, params.bootstrap_spec);
		delete alignment;
		alignment = bootstrap_alignment;
        iqtree->aln = alignment;
	}
	if (alignment->isSuperAlignment()) {
		((SuperAlignment*)alignment)->printCombinedAlignment(params.aln_output);
		if (params.print_subaln)
			((SuperAlignment*)alignment)->printSubAlignments(params, ((PhyloSuperTree*)iqtree)->part_info);

	} else if (params.gap_masked_aln) {
		Alignment out_aln;
		Alignment masked_aln(params.gap_masked_aln, params.sequence_type, params.intype);
		out_aln.createGapMaskedAlignment(&masked_aln, alignment);
		out_aln.printPhylip(params.aln_output, false, params.aln_site_list,
				params.aln_nogaps, params.aln_no_const_sites, params.ref_seq_name);
		string str = params.gap_masked_aln;
		str += ".sitegaps";
		out_aln.printSiteGaps(str.c_str());
	} else if (params.aln_output_format == ALN_PHYLIP)
		alignment->printPhylip(params.aln_output, false, params.aln_site_list,
				params.aln_nogaps, params.aln_no_const_sites, params.ref_seq_name);
	else if (params.aln_output_format == ALN_FASTA)
		alignment->printFasta(params.aln_output, false, params.aln_site_list,
				params.aln_nogaps, params.aln_no_const_sites, params.ref_seq_name);
}

/**
    2016-08-04: compute a site frequency model for profile mixture model
*/
void computeSiteFrequencyModel(Params &params, Alignment *alignment) {

    cout << endl << "===> COMPUTING SITE FREQUENCY MODEL BASED ON TREE FILE " << params.tree_freq_file << endl;
    assert(params.tree_freq_file);
    PhyloTree *tree = new PhyloTree(alignment);
    tree->setParams(&params);
    bool myrooted = params.is_rooted;
    tree->readTree(params.tree_freq_file, myrooted);
    tree->setAlignment(alignment);
    tree->setRootNode(params.root);
    
	ModelsBlock *models_block = readModelsDefinition(params);
    tree->setModelFactory(new ModelFactory(params, tree, models_block));
    delete models_block;
    tree->setModel(tree->getModelFactory()->model);
    tree->setRate(tree->getModelFactory()->site_rate);
    tree->setLikelihoodKernel(params.SSE);

    if (!tree->getModel()->isMixture())
        outError("No mixture model was specified!");
    uint64_t mem_size = tree->getMemoryRequired();
    uint64_t total_mem = getMemorySize();
    cout << "NOTE: " << (mem_size / 1024) / 1024 << " MB RAM is required!" << endl;
    if (mem_size >= total_mem) {
        outError("Memory required exceeds your computer RAM size!");
    }
#ifdef BINARY32
    if (mem_size >= 2000000000) {
        outError("Memory required exceeds 2GB limit of 32-bit executable");
    }
#endif

    tree->initializeAllPartialLh();
    tree->getModelFactory()->optimizeParameters(params.fixed_branch_length, true, params.modeps);

    size_t nptn = alignment->getNPattern(), nstates = alignment->num_states;
    double *ptn_state_freq = new double[nptn*nstates];
    tree->computePatternStateFreq(ptn_state_freq);
    alignment->site_state_freq.resize(nptn);
    for (size_t ptn = 0; ptn < nptn; ptn++) {
        double *f = new double[nstates];
        memcpy(f, ptn_state_freq+ptn*nstates, sizeof(double)*nstates);
        alignment->site_state_freq[ptn] = f;
    }
    alignment->getSitePatternIndex(alignment->site_model);
    printSiteStateFreq(((string)params.out_prefix+".sitefreq").c_str(), tree, ptn_state_freq);
    params.print_site_state_freq = WSF_NONE;
    
    delete [] ptn_state_freq;
    delete tree;
    
    cout << endl << "===> CONTINUE ANALYSIS USING THE INFERRED SITE FREQUENCY MODEL" << endl;
}


/**********************************************************
 * TOP-LEVEL FUNCTION
 ***********************************************************/
void runPhyloAnalysis(Params &params, Checkpoint *checkpoint) {
	Alignment *alignment;
	IQTree *tree;
    
    checkpoint->putBool("finished", false);
    checkpoint->setDumpInterval(params.checkpoint_dump_interval);

	/****************** read in alignment **********************/
	if (params.partition_file) {
		// Partition model analysis
		if(params.partition_type){
			// since nni5 does not work yet, stop the programm
/*			if(params.nni5)
				outError("-nni5 option is unsupported yet for proportitional partition model. please use -nni1 option");*/
//			if(params.aLRT_replicates || params.localbp_replicates)
//				outError("-alrt or -lbp option is unsupported yet for joint/proportional partition model");
			// initialize supertree - Proportional Edges case, "-spt p" option
			tree = new PhyloSuperTreePlen(params);
		} else {
			// initialize supertree stuff if user specifies partition file with -sp option
			tree = new PhyloSuperTree(params);
		}
		// this alignment will actually be of type SuperAlignment
		alignment = tree->aln;
	} else {
		alignment = new Alignment(params.aln_file, params.sequence_type, params.intype);
		if (params.freq_const_patterns) {
			int orig_nsite = alignment->getNSite();
			alignment->addConstPatterns(params.freq_const_patterns);
			cout << "INFO: " << alignment->getNSite() - orig_nsite << " const sites added into alignment" << endl;
		}
        
        if (params.tree_freq_file) {
            if (checkpoint->getBool("finishedSiteFreqFile")) {
                alignment->readSiteStateFreq(((string)params.out_prefix + ".sitefreq").c_str());
                params.print_site_state_freq = WSF_NONE;
                cout << "CHECKPOINT: Site frequency model restored" << endl;
            } else {
                computeSiteFrequencyModel(params, alignment);
                checkpoint->putBool("finishedSiteFreqFile", true);
                checkpoint->dump();
            }
        }
        if (params.site_freq_file) {
            alignment->readSiteStateFreq(params.site_freq_file);
        }
            
		tree = new IQTree(alignment);
	}

    tree->setCheckpoint(checkpoint);
    if (params.min_branch_length <= 0.0) {
        params.min_branch_length = 1e-6;
        if (tree->getAlnNSite() >= 100000) {
            params.min_branch_length = 0.1 / (tree->getAlnNSite());
            tree->num_precision = max((int)ceil(-log10(Params::getInstance().min_branch_length))+1, 6);
            cout.precision(12);
            cout << "NOTE: minimal branch length is reduced to " << params.min_branch_length << " for long alignment" << endl;
            cout.precision(3);
        }
    }

	string original_model = params.model_name;

	if (params.concatenate_aln) {
		Alignment aln(params.concatenate_aln, params.sequence_type, params.intype);
		cout << "Concatenating " << params.aln_file << " with " << params.concatenate_aln << " ..." << endl;
		alignment->concatenateAlignment(&aln);
	}

    if (params.constraint_tree_file) {
        cout << "Reading constraint tree " << params.constraint_tree_file << "..." << endl;
        tree->constraintTree.initConstraint(params.constraint_tree_file, alignment->getSeqNames());
        if (params.start_tree == STT_PLL_PARSIMONY)
            params.start_tree = STT_PARSIMONY;
        else if (params.start_tree == STT_BIONJ)
            outError("Constraint tree does not work with -t BIONJ");
            
    }

    if (params.compute_seq_identity_along_tree) {
        if (!params.user_file)
            outError("Please supply a user tree file!");
        tree->readTree(params.user_file, params.is_rooted);
        if (!tree->rooted && !params.root) {
            outError("Tree is unrooted, thus you have to specify a root with -o option");
        }
        tree->setAlignment(tree->aln);
        if (!tree->rooted)
            tree->setRootNode(params.root);
        tree->computeSeqIdentityAlongTree();
        if (verbose_mode >= VB_MED)
            tree->drawTree(cout);
        string out_tree = (string)params.out_prefix + ".seqident_tree";
        tree->printTree(out_tree.c_str());
        cout << "Tree with sequence identity printed to " << out_tree << endl;
	} else if (params.aln_output) {
		/************ convert alignment to other format and write to output file *************/
		convertAlignment(params, tree);
	} else if (params.gbo_replicates > 0 && params.user_file && params.second_tree) {
		// run one of the UFBoot analysis
//		runGuidedBootstrap(params, alignment, *tree);
		outError("Obsolete feature");
	} else if (params.avh_test) {
		// run one of the wondering test for Arndt
//		runAvHTest(params, alignment, *tree);
		outError("Obsolete feature");
	} else if (params.bootlh_test) {
		// run Arndt's plot of tree likelihoods against bootstrap alignments
//		runBootLhTest(params, alignment, *tree);
		outError("Obsolete feature");
	} else if (params.num_bootstrap_samples == 0) {
		// the main Maximum likelihood tree reconstruction
		vector<ModelInfo> *model_info = new vector<ModelInfo>;
		alignment->checkGappySeq(params.remove_empty_seq);

		// remove identical sequences
        if (params.ignore_identical_seqs) {
            tree->removeIdenticalSeqs(params);
            if (tree->removed_seqs.size() > 0) {
                string filename = (string)params.out_prefix + ".uniqueseq.phy";
                if (tree->isSuperTree())
                    ((SuperAlignment*)tree->aln)->printCombinedAlignment(filename.c_str());
                else
                    tree->aln->printPhylip(filename.c_str());
                cout << endl << "For your convenience alignment with unique sequences printed to " << filename << endl;
            }
        }
        alignment = NULL; // from now on use tree->aln instead

		// call main tree reconstruction
        runTreeReconstruction(params, original_model, *tree, *model_info);
		if (params.gbo_replicates && params.online_bootstrap) {
			if (params.print_ufboot_trees)
				tree->writeUFBootTrees(params);

			cout << endl << "Computing bootstrap consensus tree..." << endl;
			string splitsfile = params.out_prefix;
			splitsfile += ".splits.nex";
			computeConsensusTree(splitsfile.c_str(), 0, 1e6, params.split_threshold,
					params.split_weight_threshold, NULL, params.out_prefix, NULL, &params);
			// now optimize branch lengths of the consensus tree
			string current_tree = tree->getTreeString();
			splitsfile = params.out_prefix;
			splitsfile += ".contree";
			tree->readTreeFile(splitsfile);
			// bug fix
//			if ((tree->sse == LK_EIGEN || tree->sse == LK_EIGEN_SSE) && !tree->isBifurcating()) {
//				cout << "NOTE: Changing to old kernel as consensus tree is multifurcating" << endl;
//                if (tree->sse == LK_EIGEN)
//                    tree->changeLikelihoodKernel(LK_NORMAL);
//                else
//                    tree->changeLikelihoodKernel(LK_SSE);
//			}

			tree->initializeAllPartialLh();
			tree->fixNegativeBranch(true);

			tree->boot_consense_logl = tree->optimizeAllBranches();
			cout << "Log-likelihood of consensus tree: " << tree->boot_consense_logl << endl;
		    tree->setRootNode(params.root);
		    tree->insertTaxa(tree->removed_seqs, tree->twin_seqs);
			tree->printTree(splitsfile.c_str(), WT_BR_LEN | WT_BR_LEN_FIXED_WIDTH | WT_SORT_TAXA | WT_NEWLINE);
			// revert the best tree
			tree->readTreeString(current_tree);
		}
		// reinsert identical sequences
		if (tree->removed_seqs.size() > 0) {
			// BUG FIX: dont use reinsertIdenticalSeqs anymore
			tree->insertTaxa(tree->removed_seqs, tree->twin_seqs);
			tree->printResultTree();
		}
		reportPhyloAnalysis(params, original_model, *tree, *model_info);
        delete model_info;
	} else {
		// the classical non-parameter bootstrap (SBS)
		if (params.model_name.find("LINK") != string::npos || params.model_name.find("MERGE") != string::npos)
			outError("-m TESTMERGE is not allowed when doing standard bootstrap. Please first\nfind partition scheme on the original alignment and use it for bootstrap analysis");
        if (alignment->getNSeq() < 4)
            outError("It makes no sense to perform bootstrap with less than 4 sequences.");
		runStandardBootstrap(params, original_model, alignment, tree);
	}

//	if (params.upper_bound) {
//			UpperBounds(&params, alignment, tree);
//	}

	if(verbose_mode >= VB_MED){
		if(tree->isSuperTree() && params.partition_type){
			((PhyloSuperTreePlen*) tree)->printNNIcasesNUM();
		}
	}
    // 2015-09-22: bug fix, move this line to before deleting tree
    alignment = tree->aln;
	delete tree;
	// BUG FIX: alignment can be changed, should delete tree->aln instead
    // 2015-09-22: THIS IS STUPID: after deleting tree, one cannot access tree->aln anymore
//	alignment = tree->aln;
	delete alignment;

    checkpoint->putBool("finished", true);
    checkpoint->dump(true);
}

void assignBranchSupportNew(Params &params) {
	if (!params.user_file)
		outError("No trees file provided");
	if (!params.second_tree)
		outError("No target tree file provided");
	cout << "Reading tree " << params.second_tree << " ..." << endl;
	MTree tree(params.second_tree, params.is_rooted);
	cout << tree.leafNum << " taxa and " << tree.branchNum << " branches" << endl;
	tree.assignBranchSupport(params.user_file);
	string str = params.second_tree;
	str += ".suptree";
	tree.printTree(str.c_str());
	cout << "Tree with assigned branch supports written to " << str << endl;
	if (verbose_mode >= VB_DEBUG)
		tree.drawTree(cout);
}



/**
 * assign split occurence frequencies from a set of input trees onto a target tree
 * NOTE: input trees must have the same taxon set
 * @param input_trees file containing NEWICK tree strings
 * @param burnin number of beginning trees to discard
 * @param max_count max number of trees to read in
 * @param target_tree the target tree
 * @param rooted TRUE if trees are rooted, false for unrooted trees
 * @param output_file file name to write output tree with assigned support values
 * @param out_prefix prefix of output file
 * @param mytree (OUT) resulting tree with support values assigned from target_tree
 * @param tree_weight_file file containing INTEGER weights of input trees
 * @param params program parameters
 */
void assignBootstrapSupport(const char *input_trees, int burnin, int max_count,
		const char *target_tree, bool rooted, const char *output_tree,
		const char *out_prefix, MExtTree &mytree, const char* tree_weight_file,
		Params *params) {
	//bool rooted = false;
	// read the tree file
	cout << "Reading tree " << target_tree << " ..." << endl;
	mytree.init(target_tree, rooted);
	// reindex the taxa in the tree to aphabetical names
	NodeVector taxa;
	mytree.getTaxa(taxa);
	sort(taxa.begin(), taxa.end(), nodenamecmp);
	int i = 0;
	for (NodeVector::iterator it = taxa.begin(); it != taxa.end(); it++) {
		(*it)->id = i++;
	}

	/*
	 string filename = params.boot_trees;
	 filename += ".nolen";
	 boot_trees.printTrees(filename.c_str(), false);
	 return;
	 */
	SplitGraph sg;
	SplitIntMap hash_ss;
	// make the taxa name
	vector<string> taxname;
	taxname.resize(mytree.leafNum);
	mytree.getTaxaName(taxname);

	// read the bootstrap tree file
	double scale = 100.0;
	if (params->scaling_factor > 0)
		scale = params->scaling_factor;

	MTreeSet boot_trees;
	if (params && detectInputFile((char*) input_trees) == IN_NEXUS) {
		sg.init(*params);
		for (SplitGraph::iterator it = sg.begin(); it != sg.end(); it++)
			hash_ss.insertSplit((*it), (*it)->getWeight());
		StrVector sgtaxname;
		sg.getTaxaName(sgtaxname);
		i = 0;
		for (StrVector::iterator sit = sgtaxname.begin();
				sit != sgtaxname.end(); sit++, i++) {
			Node *leaf = mytree.findLeafName(*sit);
			if (!leaf)
				outError("Tree does not contain taxon ", *sit);
			leaf->id = i;
		}
		scale /= sg.maxWeight();
	} else {
		boot_trees.init(input_trees, rooted, burnin, max_count,
				tree_weight_file);
		boot_trees.convertSplits(taxname, sg, hash_ss, SW_COUNT, -1, params->support_tag);
		scale /= boot_trees.sumTreeWeights();
	}
	//sg.report(cout);
	cout << "Rescaling split weights by " << scale << endl;
	if (params->scaling_factor < 0)
		sg.scaleWeight(scale, true);
	else {
		sg.scaleWeight(scale, false, params->numeric_precision);
	}

	cout << sg.size() << " splits found" << endl;
	// compute the percentage of appearance
	//	printSplitSet(sg, hash_ss);
	//sg.report(cout);
	cout << "Creating bootstrap support values..." << endl;
	mytree.createBootstrapSupport(taxname, boot_trees, sg, hash_ss, params->support_tag);
	//mytree.scaleLength(100.0/boot_trees.size(), true);
	string out_file;
	if (output_tree)
		out_file = output_tree;
	else {
		if (out_prefix)
			out_file = out_prefix;
		else
			out_file = target_tree;
		out_file += ".suptree";
	}

	mytree.printTree(out_file.c_str());
	cout << "Tree with assigned bootstrap support written to " << out_file
			<< endl;
	/*
	if (out_prefix)
		out_file = out_prefix;
	else
		out_file = target_tree;
	out_file += ".supval";
	mytree.writeInternalNodeNames(out_file);

	cout << "Support values written to " << out_file << endl;
	*/
}

void computeConsensusTree(const char *input_trees, int burnin, int max_count,
		double cutoff, double weight_threshold, const char *output_tree,
		const char *out_prefix, const char *tree_weight_file, Params *params) {
	bool rooted = false;

	// read the bootstrap tree file
	/*
	 MTreeSet boot_trees(input_trees, rooted, burnin, tree_weight_file);
	 string first_taxname = boot_trees.front()->root->name;
	 //if (params.root) first_taxname = params.root;

	 SplitGraph sg;

	 boot_trees.convertSplits(sg, cutoff, SW_COUNT, weight_threshold);*/

	//sg.report(cout);
	SplitGraph sg;
	SplitIntMap hash_ss;
	// make the taxa name
	//vector<string> taxname;
	//taxname.resize(mytree.leafNum);
	//mytree.getTaxaName(taxname);

	// read the bootstrap tree file
	double scale = 100.0;
	if (params->scaling_factor > 0)
		scale = params->scaling_factor;

	MTreeSet boot_trees;
	if (params && detectInputFile((char*) input_trees) == IN_NEXUS) {
		char *user_file = params->user_file;
		params->user_file = (char*) input_trees;
		params->split_weight_summary = SW_COUNT; // count number of splits
		sg.init(*params);
		params->user_file = user_file;
		for (SplitGraph::iterator it = sg.begin(); it != sg.end(); it++)
			hash_ss.insertSplit((*it), (*it)->getWeight());
		/*		StrVector sgtaxname;
		 sg.getTaxaName(sgtaxname);
		 i = 0;
		 for (StrVector::iterator sit = sgtaxname.begin(); sit != sgtaxname.end(); sit++, i++) {
		 Node *leaf = mytree.findLeafName(*sit);
		 if (!leaf) outError("Tree does not contain taxon ", *sit);
		 leaf->id = i;
		 }*/
		scale /= sg.maxWeight();
	} else {
		boot_trees.init(input_trees, rooted, burnin, max_count,
				tree_weight_file);
		boot_trees.convertSplits(sg, cutoff, SW_COUNT, weight_threshold);
		scale /= boot_trees.sumTreeWeights();
		cout << sg.size() << " splits found" << endl;
	}
	//sg.report(cout);
	if (verbose_mode >= VB_MED)
		cout << "Rescaling split weights by " << scale << endl;
	if (params->scaling_factor < 0)
		sg.scaleWeight(scale, true);
	else {
		sg.scaleWeight(scale, false, params->numeric_precision);
	}



	//cout << "Creating greedy consensus tree..." << endl;
	MTree mytree;
	SplitGraph maxsg;
	sg.findMaxCompatibleSplits(maxsg);

	if (verbose_mode >= VB_MAX)
		maxsg.saveFileStarDot(cout);
	//cout << "convert compatible split system into tree..." << endl;
	mytree.convertToTree(maxsg);
	//cout << "done" << endl;
	string taxname;
	if (params->root)
		taxname = params->root;
	else
		taxname = sg.getTaxa()->GetTaxonLabel(0);
	Node *node = mytree.findLeafName(taxname);
	if (node)
		mytree.root = node;
	// mytree.scaleLength(100.0 / boot_trees.sumTreeWeights(), true);

	// mytree.getTaxaID(maxsg.getSplitsBlock()->getCycle());
	//maxsg.saveFile(cout);

	string out_file;

	if (output_tree)
		out_file = output_tree;
	else {
		if (out_prefix)
			out_file = out_prefix;
		else
			out_file = input_trees;
		out_file += ".contree";
	}

//	if (removed_seqs.size() > 0)
//		mytree.insertTaxa(removed_seqs, twin_seqs);

	mytree.printTree(out_file.c_str(), WT_BR_CLADE);
	cout << "Consensus tree written to " << out_file << endl;

	if (output_tree)
		out_file = output_tree;
	else {
		if (out_prefix)
			out_file = out_prefix;
		else
			out_file = input_trees;
		out_file += ".splits";
	}

    //sg.scaleWeight(0.01, false, 4);
	if (params->print_splits_file) {
		sg.saveFile(out_file.c_str(), IN_OTHER, true);
		cout << "Non-trivial split supports printed to star-dot file " << out_file << endl;
	}

}

void computeConsensusNetwork(const char *input_trees, int burnin, int max_count,
		double cutoff, int weight_summary, double weight_threshold, const char *output_tree,
		const char *out_prefix, const char* tree_weight_file) {
	bool rooted = false;

	// read the bootstrap tree file
	MTreeSet boot_trees(input_trees, rooted, burnin, max_count,
			tree_weight_file);

	SplitGraph sg;
	//SplitIntMap hash_ss;

	boot_trees.convertSplits(sg, cutoff, weight_summary, weight_threshold);

	string out_file;

	if (output_tree)
		out_file = output_tree;
	else {
		if (out_prefix)
			out_file = out_prefix;
		else
			out_file = input_trees;
		out_file += ".nex";
	}

	sg.saveFile(out_file.c_str(), IN_NEXUS);
	cout << "Consensus network printed to " << out_file << endl;

	if (output_tree)
		out_file = output_tree;
	else {
		if (out_prefix)
			out_file = out_prefix;
		else
			out_file = input_trees;
		out_file += ".splits";
	}
	if (verbose_mode >= VB_MED) {
		sg.saveFile(out_file.c_str(), IN_OTHER, true);
		cout << "Non-trivial split supports printed to star-dot file " << out_file << endl;
	}

}

