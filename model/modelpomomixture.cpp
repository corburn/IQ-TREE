//
//  modelpomomixture.cpp
//  iqtree
//
//  Created by Minh Bui on 7/22/16.
//
//

#include "modelpomomixture.h"
#include "rategamma.h"

ModelPoMoMixture::ModelPoMoMixture(const char *model_name,
                     string model_params,
                     StateFreqType freq_type,
                     string freq_params,
                     PhyloTree *tree,
                     bool is_reversible,
                     string pomo_params, string pomo_rate_str)
	:  
        ModelGTR(tree),
        ModelPoMo(model_name, model_params, freq_type, freq_params, tree, is_reversible, pomo_params),
        ModelMixture(tree)
    
{
    // get number of categories
    int m, num_rate_cats = 4;
    if (pomo_rate_str.length() > 2 && isdigit(pomo_rate_str[2])) {
        int end_pos;
        num_rate_cats = convert_int(pomo_rate_str.substr(2).c_str(), end_pos);
        if (num_rate_cats < 1) outError("Wrong number of rate categories");
    }

    // initialize rate heterogeneity
    ratehet = new RateGamma(num_rate_cats, Params::getInstance().gamma_shape, Params::getInstance().gamma_median, tree);

    // initialize mixture
    prop = aligned_alloc<double>(num_rate_cats);
    
    for (m = 0; m < num_rate_cats; m++) {
        ModelGTR* model = new ModelGTR(tree);
        model->init(FREQ_USER_DEFINED);
        model->total_num_subst = ratehet->getRate(m);
        push_back(model);
        prop[m] = ratehet->getProp(m);
    }
    initMem();


    ModelGTR::init(FREQ_USER_DEFINED);
    
}




ModelPoMoMixture::~ModelPoMoMixture() {
}

void ModelPoMoMixture::saveCheckpoint() {
    ModelPoMo::saveCheckpoint();
}

void ModelPoMoMixture::restoreCheckpoint() {
    ModelPoMo::restoreCheckpoint();
}


int ModelPoMoMixture::getNDim() {
    return ModelPoMo::getNDim();
}


int ModelPoMoMixture::getNDimFreq() {
    return ModelPoMo::getNDimFreq();
}


double ModelPoMoMixture::targetFunk(double x[]) {
    return ModelPoMo::targetFunk(x);
}



void ModelPoMoMixture::setBounds(double *lower_bound, double *upper_bound, bool *bound_check) {
    ModelPoMo::setBounds(lower_bound, upper_bound, bound_check);
}


void ModelPoMoMixture::writeInfo(ostream &out) {
    ModelPoMo::writeInfo(out);
}

void ModelPoMoMixture::decomposeRateMatrix() {
    ModelPoMo::decomposeRateMatrix();
    // propagate eigenvalues and eigenvectors
    int m, nmix = getNMixtures(), num_states_2 = num_states*num_states;
    for (m = 1; m < nmix; m++) {
        memcpy(eigenvalues+m*num_states, eigenvalues+(m-1)*num_states, sizeof(double)*num_states);
        memcpy(eigenvectors+m*num_states_2, eigenvectors+(m-1)*num_states_2, sizeof(double)*num_states_2);
        memcpy(inv_eigenvectors+m*num_states_2, inv_eigenvectors+(m-1)*num_states_2, sizeof(double)*num_states_2);
    }
}


void ModelPoMoMixture::setVariables(double *variables) {
    ModelPoMo::setVariables(variables);
}


bool ModelPoMoMixture::getVariables(double *variables) {
    return ModelPoMo::getVariables(variables);
}


double ModelPoMoMixture::optimizeParameters(double gradient_epsilon) {
    return ModelPoMo::optimizeParameters(gradient_epsilon);
}