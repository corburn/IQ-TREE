/*
 * phylotreeavx.cpp
 *
 *  Created on: Dec 14, 2014
 *      Author: minh
 */


#include "phylokernel.h"
#include "phylokernelmixture.h"
#include "phylokernelmixrate.h"
#include "vectorclass/vectorclass.h"

#ifndef __AVX__
#error "You must compile this file with AVX enabled!"
#endif

void PhyloTree::setDotProductAVX() {
#ifdef BOOT_VAL_FLOAT
		dotProduct = &PhyloTree::dotProductSIMD<float, Vec8f, 8>;
#else
		dotProduct = &PhyloTree::dotProductSIMD<double, Vec4d, 4>;
#endif

}

void PhyloTree::setLikelihoodKernelAVX() {
	switch(aln->num_states) {
	case 4:
		if (model_factory && model_factory->model->isMixture()) {
			if (model_factory->fused_mix_rate) {
				computeLikelihoodBranchPointer = &PhyloTree::computeMixrateLikelihoodBranchEigenSIMD<Vec4d, 4, 4>;
				computeLikelihoodDervPointer = &PhyloTree::computeMixrateLikelihoodDervEigenSIMD<Vec4d, 4, 4>;
				computePartialLikelihoodPointer = &PhyloTree::computeMixratePartialLikelihoodEigenSIMD<Vec4d, 4, 4>;
				computeLikelihoodFromBufferPointer = &PhyloTree::computeMixrateLikelihoodFromBufferEigenSIMD<Vec4d, 4, 4>;
			} else {
				computeLikelihoodBranchPointer = &PhyloTree::computeMixtureLikelihoodBranchEigenSIMD<Vec4d, 4, 4>;
				computeLikelihoodDervPointer = &PhyloTree::computeMixtureLikelihoodDervEigenSIMD<Vec4d, 4, 4>;
				computePartialLikelihoodPointer = &PhyloTree::computeMixturePartialLikelihoodEigenSIMD<Vec4d, 4, 4>;
				computeLikelihoodFromBufferPointer = &PhyloTree::computeMixtureLikelihoodFromBufferEigenSIMD<Vec4d, 4, 4>;
			}
		} else {
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec4d, 4, 4>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec4d, 4, 4>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec4d, 4, 4>;
			computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 4>;
		}
		break;
	case 20:
		if (model_factory && model_factory->model->isMixture()) {
			if (model_factory->fused_mix_rate) {
				computeLikelihoodBranchPointer = &PhyloTree::computeMixrateLikelihoodBranchEigenSIMD<Vec4d, 4, 20>;
				computeLikelihoodDervPointer = &PhyloTree::computeMixrateLikelihoodDervEigenSIMD<Vec4d, 4, 20>;
				computePartialLikelihoodPointer = &PhyloTree::computeMixratePartialLikelihoodEigenSIMD<Vec4d, 4, 20>;
				computeLikelihoodFromBufferPointer = &PhyloTree::computeMixrateLikelihoodFromBufferEigenSIMD<Vec4d, 4, 20>;
			} else {
				computeLikelihoodBranchPointer = &PhyloTree::computeMixtureLikelihoodBranchEigenSIMD<Vec4d, 4, 20>;
				computeLikelihoodDervPointer = &PhyloTree::computeMixtureLikelihoodDervEigenSIMD<Vec4d, 4, 20>;
				computePartialLikelihoodPointer = &PhyloTree::computeMixturePartialLikelihoodEigenSIMD<Vec4d, 4, 20>;
				computeLikelihoodFromBufferPointer = &PhyloTree::computeMixtureLikelihoodFromBufferEigenSIMD<Vec4d, 4, 20>;
			}
		} else {
			computeLikelihoodBranchPointer = &PhyloTree::computeLikelihoodBranchEigenSIMD<Vec4d, 4, 20>;
			computeLikelihoodDervPointer = &PhyloTree::computeLikelihoodDervEigenSIMD<Vec4d, 4, 20>;
			computePartialLikelihoodPointer = &PhyloTree::computePartialLikelihoodEigenSIMD<Vec4d, 4, 20>;
			computeLikelihoodFromBufferPointer = &PhyloTree::computeLikelihoodFromBufferEigenSIMD<Vec4d, 4, 20>;
		}
		break;
	default:
		assert(0);
		break;
	}
}
