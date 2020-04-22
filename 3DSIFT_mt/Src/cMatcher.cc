
#include"../../Include/cMatcher.h"

#include<iostream>
#include<sstream>
#include<fstream>
#include<time.h>
#include<cfloat>
#include<Eigen/core>
#include<Eigen/dense>

using namespace std;

namespace CPUSIFT{

	//utils funcs
	double KP_squareSum(const Keypoint &k1, const Keypoint &k2) {
		double result = 0;
		for (int i = 0; i < DESC_LENGTH; i++) {
			result += k1.desc[i] * k2.desc[i];
		}
		return result;
	}

	muBruteMatcher::muBruteMatcher() {
		thread_num = omp_get_num_procs();
		cout << "matching using num of threads:" << thread_num << endl;
	}

	float muBruteMatcher::getCalculationTime() { return totalTime; };

	std::vector<float> muBruteMatcher::getGlodenDistSquare() { return glodenDistSquare; };

	std::vector<float> muBruteMatcher::getSilverDistSquare() { return silverDistSquare; };

	std::vector<int> muBruteMatcher::getGlodenIdx() { return glodenIdx; };

	std::vector<int> muBruteMatcher::getSilverIdx() { return silverIdx; };

	void muBruteMatcher::calMatches(vector<float> &gDistSquare_, vector<float> &sDistSquare_, vector<int> &gIdx_, vector<int> &sIdx_, const vector<Keypoint>& ref_kp_, const vector<Keypoint>& tar_kp_, const vector<int> *mask) {

		const int num1 = ref_kp_.size();
		const int num2 = tar_kp_.size();
		//calculation
#pragma omp parallel for schedule(dynamic) num_threads(thread_num) 
		for (int i = 0; i < num1; i++)
		{
			if (mask != nullptr && (*mask)[i] == 0) {
				//masked assigned with -1, representing invalid matches
				gIdx_[i] = -1;
				continue;
			}
			auto &kpr = ref_kp_[i];
			double d1 = FLT_MIN;
			double d2 = FLT_MIN;
			int i1 = -1, i2 = -1;

			for (int j = 0; j < num2; j++) {
				double sij = KP_squareSum(kpr, tar_kp_[j]);
				if (sij > d1) {
					d2 = d1;
					i2 = i1;
					d1 = sij;
					i1 = j;
				}
				else if (sij > d2) {
					d2 = sij;
					i2 = j;
				}
			}
			d2 = 2 - 2 * d2;
			d1 = 2 - 2 * d1;

			gDistSquare_[i] = d1;
			sDistSquare_[i] = d2;
			gIdx_[i] = i1;
			sIdx_[i] = i2;
		}
	}

	void muBruteMatcher::filter(vector<int> &gIdx_, const vector<float> &gDistSquare_, const vector<float> &sDistSquare_, const double thresHold) {
		const int num_ = gDistSquare_.size();
		const double thresSquare = thresHold * thresHold;

#pragma omp parallel for num_threads(thread_num) 
		for (int i = 0; i < num_; ++i) {
			//only valid matches are judged
			if (gIdx_[i] < 0)
				continue;
			auto d1 = gDistSquare_[i];
			auto d2 = sDistSquare_[i];
			if (d1 / d2 >= thresSquare) {
				gIdx_[i] *= -1;
			}
		}
		return;
	}

	void muBruteMatcher::toCvec(vector<Cvec>& refMatch, vector<Cvec>& tarMatch, const vector<Keypoint>& ref_kp, const vector<Keypoint>& tar_kp, const vector<int> &gIdx_) {

		const int num = ref_kp.size();
		for (int i = 0; i < num; ++i) {
			int j = gIdx_[i];
			if (j < 0)
				continue;
			Cvec ref_ = { ref_kp[i].rx, ref_kp[i].ry, ref_kp[i].rz };
			Cvec tar_ = { tar_kp[j].rx, tar_kp[j].ry, tar_kp[j].rz };
			refMatch.push_back(ref_);
			tarMatch.push_back(tar_);
		}
		return;
	}

	void muBruteMatcher::countMatched(vector<int> &countMatched_, const vector<int> &gIdx_) {
		for (int i = 0; i < gIdx_.size(); ++i) {
			int idx = gIdx_[i];
			if (idx >= 0)
				countMatched_[idx] += 1;
		}
	}

	vector<int> muBruteMatcher::toMask(const vector<int> &countMatched_, const int countThres_) {
		auto mask = countMatched_;
		for (auto &i : mask) {
			if (i > countThres_)
				i = 1;
			else
				i = 0;
		}
		return mask;
	}

	void muBruteMatcher::bijectFilter(std::vector<int> &refIdx_, const std::vector<int> &tarMask_, const std::vector<int> &tarIdx_) {
		const int num = refIdx_.size();
#pragma omp parallel for num_threads(thread_num) 
		for (int i = 0; i < num; ++i) {
			int mIdx = refIdx_[i];
			if (mIdx < 0 || tarMask_[mIdx] == 0)
				continue;
			//reject non-biject 
			if (tarIdx_[mIdx] != i)
				refIdx_[i] *= -1;
		}
	}

	void muBruteMatcher::bijectMatchBase(vector<Cvec>& refMatch, vector<Cvec>& tarMatch, const vector<Keypoint>& ref_kp, const vector<Keypoint>& tar_kp, const double thresHold, const matchType type) {

		auto tStart = omp_get_wtime();
		const int num1 = ref_kp.size();
		const int num2 = tar_kp.size();

		glodenDistSquare.resize(num1);
		silverDistSquare.resize(num1);
		glodenIdx = vector<int>(num1, -1);
		silverIdx = vector<int>(num1, -1);

		glodenDistSquare2.resize(num2);
		silverDistSquare2.resize(num2);
		glodenIdx2 = vector<int>(num2, -1);
		silverIdx2 = vector<int>(num2, -1);
		matchedCounts = vector<int>(num2, 0);

		double tMatchEnd = 0.0, tFilterEnd = 0.0, tCountMatchEnd = 0.0, tRevMatchEnd = 0.0, tRevFilterEnd = 0.0, tBijectFilterEnd = 0.0, tConverseEnd = 0.0;

		//match
		//-----------------------Ref -> Tar-----------------------
		calMatches(glodenDistSquare, silverDistSquare, glodenIdx, silverIdx, ref_kp, tar_kp);
		tMatchEnd = omp_get_wtime();
		//filter
		filter(glodenIdx, glodenDistSquare, silverDistSquare, thresHold);
		tFilterEnd = omp_get_wtime();


		if (type != inject) {
			const int maskThres = (type == biject ? 0 : 1);
			//-----------------------Ref <- Tar-----------------------
			//count match, the maskThres is chosen according to matchType
			countMatched(matchedCounts, glodenIdx);
			vector<int> mask = toMask(matchedCounts, maskThres);
			tCountMatchEnd = omp_get_wtime();
			//rev match
			calMatches(glodenDistSquare2, silverDistSquare2, glodenIdx2, silverIdx2, tar_kp, ref_kp, &mask);
			tRevMatchEnd = omp_get_wtime();

			//rev filter if biject
			//if (type == biject) {
				filter(glodenIdx2, glodenDistSquare2, silverDistSquare2, thresHold);
			//}
			tRevFilterEnd = omp_get_wtime();

			//biject filter
			//-------------------Ref <-match-> Tar-------------------
			bijectFilter(glodenIdx, mask, glodenIdx2);
			tBijectFilterEnd = omp_get_wtime();

			//
			countMatchedTime = tCountMatchEnd - tFilterEnd;
			revMatchTime = tRevMatchEnd - tCountMatchEnd;
			revFilterTime = tRevFilterEnd - tRevMatchEnd;
			bijectFilterTime = tBijectFilterEnd - tRevFilterEnd;
		}
		else {
			countMatchedTime = revMatchTime = revFilterTime = bijectFilterTime = 0.0;
		}

		//converse
		toCvec(refMatch, tarMatch, ref_kp, tar_kp, glodenIdx);
		tConverseEnd = omp_get_wtime();

		matchTime = tMatchEnd - tStart;
		filterTime = tFilterEnd - tMatchEnd;
		converseTime = tConverseEnd - (type == inject ? tFilterEnd : tBijectFilterEnd);
		totalTime = tConverseEnd - tStart;
		return;
	}

	//public
	void muBruteMatcher::injectMatch(vector<Cvec>& refMatch, vector<Cvec>& tarMatch, const vector<Keypoint>& ref_kp, const vector<Keypoint>& tar_kp, const double thresHold) {
		bijectMatchBase(refMatch, tarMatch, ref_kp, tar_kp, thresHold, inject);
	}

	void muBruteMatcher::bijectMatch(vector<Cvec>& refMatch, vector<Cvec>& tarMatch, const vector<Keypoint>& ref_kp, const vector<Keypoint>& tar_kp, const double thresHold) {
		bijectMatchBase(refMatch, tarMatch, ref_kp, tar_kp, thresHold, biject);
	}

	void muBruteMatcher::enhancedMatch(vector<Cvec>& refMatch, vector<Cvec>& tarMatch, const vector<Keypoint>& ref_kp, const vector<Keypoint>& tar_kp, const double thresHold) {
		bijectMatchBase(refMatch, tarMatch, ref_kp, tar_kp, thresHold, enhanced);
	}

}