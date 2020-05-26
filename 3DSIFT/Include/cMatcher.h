#ifndef _CPU_MATCHER_H_
#define _CPU_MATCHER_H_

#include <vector>
#include "cSIFT3D.h"
#include "Util/common.h"

#define DESC_LENGTH 768

namespace CPUSIFT {

	class muBruteMatcher {
	private:
		enum matchType {
			inject = 1,
			biject = 2,
			enhanced = 3
		};
		int thread_num = 20;
		std::vector<float> glodenDistSquare;
		std::vector<float> silverDistSquare;
		std::vector<int> glodenIdx;
		std::vector<int> silverIdx;

		//for target features set
		std::vector<int> matchedCounts;
		std::vector<float> glodenDistSquare2;
		std::vector<float> silverDistSquare2;
		std::vector<int> glodenIdx2;
		std::vector<int> silverIdx2;

		//common utils
		void calMatches(
			std::vector<float> &gDistSquare_, std::vector<float> &sDistSquare_, std::vector<int> &gIdx_, std::vector<int> &sIdx_,
			const std::vector<CPUSIFT::Keypoint> &ref_kp_, const std::vector<CPUSIFT::Keypoint> &tar_kp_,
			const std::vector<int> *mask = nullptr);

		void filter(
			std::vector<int> &gIdx_, const std::vector<float> &gDistSquare_, const std::vector<float> &sDistSquare_,
			const double thresHold);

		void toCvec(
			std::vector<CPUSIFT::Cvec>& refMatch, std::vector<CPUSIFT::Cvec> &tarMatch,
			const std::vector<CPUSIFT::Keypoint> &ref_kp, const std::vector<CPUSIFT::Keypoint> &tar_kp, const std::vector<int> &gIdx_);

		//using for biject matching
		void countMatched(std::vector<int> &countMatched_, const std::vector<int> &gIdx_);

		std::vector<int> toMask(const std::vector<int> &countMatched_, const int countThres_);

		void bijectFilter(std::vector<int> &refIdx_, const std::vector<int> &tarMask_, const std::vector<int> &tarIdx_);

		//basic matcher
		void bijectMatchBase(
			std::vector<CPUSIFT::Cvec> &refMatch, std::vector<CPUSIFT::Cvec> &tarMatch,
			const std::vector<CPUSIFT::Keypoint> &ref_kp, const std::vector<CPUSIFT::Keypoint> &tar_kp,
			const double thresHold, const matchType type);

	public:
		float matchTime = 0.0;
		float filterTime = 0.0;
		float countMatchedTime = 0.0;
		float revMatchTime = 0.0;
		float revFilterTime = 0.0;
		float bijectFilterTime = 0.0;
		float converseTime = 0.0;
		float totalTime = 0.0;

		SIFT_LIBRARY_API muBruteMatcher();
		SIFT_LIBRARY_API float getCalculationTime();
		SIFT_LIBRARY_API std::vector<float> getGlodenDistSquare();
		SIFT_LIBRARY_API std::vector<float> getSilverDistSquare();
		SIFT_LIBRARY_API std::vector<int> getGlodenIdx();
		SIFT_LIBRARY_API std::vector<int> getSilverIdx();

		SIFT_LIBRARY_API void injectMatch(
			std::vector<CPUSIFT::Cvec> &refMatch, std::vector<CPUSIFT::Cvec> &tarMatch,
			const std::vector<CPUSIFT::Keypoint> &ref_kp, const std::vector<CPUSIFT::Keypoint> &tar_kp, const double thresHold = 0.85);

		SIFT_LIBRARY_API void bijectMatch(
			std::vector<CPUSIFT::Cvec> &refMatch, std::vector<CPUSIFT::Cvec> &tarMatch,
			const std::vector<CPUSIFT::Keypoint>& ref_kp, const std::vector<CPUSIFT::Keypoint> &tar_kp, const double thresHold = 0.85);

		SIFT_LIBRARY_API void enhancedMatch(
			std::vector<CPUSIFT::Cvec>& refMatch, std::vector<CPUSIFT::Cvec> &tarMatch,
			const std::vector<CPUSIFT::Keypoint> &ref_kp, const std::vector<CPUSIFT::Keypoint> &tar_kp, const double thresHold = 0.85);

	};

}

#endif 
