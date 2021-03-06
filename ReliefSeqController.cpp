/* 
 * ReliefSeqController.cpp - Bill White - 12/20/12
 */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <omp.h>

#include <gsl/gsl_rng.h>

#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/progress.hpp>

// class interface
#include "ReliefSeqController.h"

// reliefseq project
#include "Insilico.h"
#include "Dataset.h"
#include "Statistics.h"
#include "StringUtils.h"
#include "ReliefF.h"
#include "RReliefF.h"
#include "ReliefFSeq.h"

using namespace std;
namespace po = boost::program_options;
using namespace boost;
using namespace insilico;

bool scoresSortAsc(const pair<double, string>& p1,
        const pair<double, string>& p2) {
  return p1.first < p2.first;
}

bool scoresSortAscByName(const pair<double, string>& p1,
        const pair<double, string>& p2) {
  return p1.second < p2.second;
}

bool scoresSortDesc(const pair<double, string>& p1,
        const pair<double, string>& p2) {
  return p1.first > p2.first;
}

ReliefSeqController::ReliefSeqController(Dataset* ds, po::variables_map& vm,
        AnalysisType anaType) {
  cout << Timestamp() << "ReliefSeq controller initialization:" << endl;
  if(ds) {
    dataset = ds;
  } else {
    cerr << "ERROR: ReliefSeqController: data set is not initialized" << endl;
    exit(EXIT_FAILURE);
  }
  paramsMap = vm;
  analysisType = anaType;

  /// set the relief algorithm
  algorithmMode = paramsMap["algorithm-mode"].as<string>();
  if(algorithmMode == "relieff") {
    if(ds->HasContinuousPhenotypes()) {
      cout << Timestamp() << "Constructing Regression ReliefF..." << endl;
      reliefseqAlgorithm = new RReliefF(ds, vm);
    } else {
      cout << Timestamp() << "Constructing Standard ReliefF..." << endl;
      reliefseqAlgorithm = new ReliefF(ds, vm, anaType);
    }
  } else {
    if(algorithmMode == "reliefseq") {
      cout << Timestamp() << "Constructing ReliefSeq..." << endl;
      reliefseqAlgorithm = new ReliefFSeq(ds, vm);
    } else {
      cerr << "ERROR: unrecognized ReliefSeq algorithm mode: "
              << algorithmMode << endl;
      exit(EXIT_FAILURE);
    }
  }

  outFilesPrefix = paramsMap["out-files-prefix"].as<string>();

  // set the number of attributes to remove per iteration
  numToRemovePerIteration = 0;
  if(paramsMap.count("iter-remove-n")) {
    numToRemovePerIteration = paramsMap["iter-remove-n"].as<unsigned int>();
  }
  if(paramsMap.count("iter-remove-percent")) {
    unsigned int iterPercentToRemove =
            paramsMap["iter-remove-percent"].as<unsigned int>();
    cout << Timestamp() << "Iteratively removing " << iterPercentToRemove
            << " percent" << endl;
    numToRemovePerIteration = (unsigned int) (((double) iterPercentToRemove
            / 100.0) * dataset->NumVariables());
  }
  cout << Timestamp() << "ReliefSeq will remove " << numToRemovePerIteration
          << " attributes on first iteration" << endl;

  // set the number of target attributes
  numTargetAttributes = 0;
  if(vm.count("num-target")) {
    numTargetAttributes = vm["num-target"].as<unsigned int>();
  }
  if(numTargetAttributes == 0) {
    numTargetAttributes = ds->NumVariables();
    numToRemovePerIteration = 0;
  }
  if(numTargetAttributes > dataset->NumVariables()) {
    cerr << "--num-target must be less than or equal to the "
            << "number of attributes in the data set" << endl;
    exit(EXIT_FAILURE);
  }
  cout << Timestamp() << "ReliefSeq is removing attributes until best "
          << numTargetAttributes << " remain" << endl;

  // multicore setup
  unsigned int maxThreads = omp_get_num_procs();
  cout << Timestamp() << maxThreads << " OpenMP processors available"
          << endl;
  numThreads = maxThreads;
  cout << Timestamp() << "ReliefSeq will use " << numThreads << " threads" << endl;

} // end of constructor

ReliefSeqController::~ReliefSeqController() {
  if(reliefseqAlgorithm) {
    delete reliefseqAlgorithm;
  }
}

bool ReliefSeqController::ComputeScores() {
  unsigned int numWorkingAttributes = dataset->NumVariables();
  if(numWorkingAttributes < numTargetAttributes) {
    cerr << "ERROR: The number of attributes in the data set "
            << numWorkingAttributes
            << " is less than the number of target attributes "
            << numTargetAttributes << endl;
    return false;
  }

  unsigned int iteration = 1;
  boost::progress_timer t;
  while(numWorkingAttributes >= numTargetAttributes) {
    pair<unsigned int, unsigned int> titvCounts =
            dataset->GetAttributeTiTvCounts();
    double titvRatio = titvCounts.first;
    if(titvCounts.second) {
      titvRatio = (double) titvCounts.first / (double) titvCounts.second;
    }

    cout << Timestamp()
            << "----------------------------------------------------"
            << "-------------------------" << endl;
    cout << Timestamp() << "Reliefseq algorithm...iteration: " << iteration
            << ", working attributes: " << numWorkingAttributes
            << ", target attributes: " << numTargetAttributes << endl;
    cout << Timestamp()
            << "Ti/Tv: transitions: " << titvCounts.first
            << " transversions: " << titvCounts.second
            << " ratio: " << titvRatio
            << endl;

    // -------------------------------------------------------------------------
    cout << Timestamp() << "Running ReliefF" << endl;
    if(!RunReliefF()) {
      cerr << "ERROR: ReliefF failed. Exiting." << endl;
      return false;
    }
    cout << Timestamp() << "ReliefF finished in " << t.elapsed() << " secs"
            << endl;
    if(numWorkingAttributes == numTargetAttributes) {
      sort(scores.begin(), scores.end(), scoresSortDesc);
      return true;
    }

    // write scores for each iteration
    stringstream scoreFilename;
    scoreFilename << "reliefseq." << iteration << ".scores.dat";
    ofstream outFile;
    outFile.open(scoreFilename.str().c_str());
    if(outFile.bad()) {
      cerr << "ERROR: Could not open scores file " << scoreFilename.str()
              << "for writing" << endl;
      exit(1);
    }
    cout << Timestamp()
            << "Writing ALL reliefseq scores to [" + scoreFilename.str() + "]" << endl;
    PrintAttributeScores(outFile);
    outFile.close();

    // -------------------------------------------------------------------------
    // remove the worst attributes and iterate
    cout << Timestamp() << "Removing the worst attributes" << endl;
    unsigned int numToRemove = numToRemovePerIteration;
    numToRemoveNextIteration = numToRemove - numToRemovePerIteration;
    if(paramsMap.count("iter-remove-percent")) {
      unsigned int iterPercentToRemove = paramsMap["iter-remove-percent"].as<
              unsigned int>();
      numToRemove = (int) (((double) iterPercentToRemove / 100.0)
              * dataset->NumVariables());
      numToRemoveNextIteration = (int) (((double) iterPercentToRemove / 100.0)
              * numToRemove);
      numToRemovePerIteration = numToRemove;
    }
    if((numWorkingAttributes - numToRemove) < numTargetAttributes) {
      numToRemove = numWorkingAttributes - numTargetAttributes;
    }
    if(numToRemove < 1) {
      //      cerr << "ERROR: Number of attributes to remove is less than one." << endl;
      //      return false;
      break;
    }
    cout << Timestamp() << "Removing the worst " << numToRemove
            << " attributes" << endl;
    if(!RemoveWorstAttributes(numToRemove)) {
      cerr << "ERROR: RemoveWorstAttribute failed" << endl;
      return false;
    }
    numWorkingAttributes -= numToRemove;
    cout << Timestamp() << "Attribute removal complete in " << t.elapsed()
            << " secs" << endl;

    ++iteration;
  }

  cout << Timestamp() << "ReliefSeq ran for " << iteration << " iterations"
          << endl;

  return true;
}

bool ReliefSeqController::SetKoptParameters() {
  // TODO: wrap these conversions in try/catch exception handler
  unsigned int tempKoptBegin = paramsMap["kopt-begin"].as<unsigned int>();
  unsigned int tempKoptEnd = paramsMap["kopt-end"].as<unsigned int>();
  unsigned int tempKoptStep = paramsMap["kopt-step"].as<unsigned int>();
  // changed for continuous phenos - 7/26/15
  unsigned int kmax = dataset->NumInstances();

  // error conditions
  if(tempKoptBegin > tempKoptEnd) {
    cerr << "ERROR: k optimization begin [" << tempKoptBegin
            << "] is greater than end [" << tempKoptEnd << "]" << endl;
    return false;
  }
  if(tempKoptEnd > kmax) {
    cerr << "ERROR: k optimization end [" << tempKoptEnd
            << "] is greater than maximum k [" << kmax << "]" << endl;
    return false;
  }
  if((tempKoptBegin == tempKoptEnd) == tempKoptStep) {
    cerr << "ERROR: k optimization specified but the range "
            << "and step values do not specify any iterations"
            << endl;
    return false;
  }

  // passed all error checks
  koptBegin = tempKoptBegin;
  koptEnd = tempKoptEnd;
  koptStep = tempKoptStep;
  cout << Timestamp() << "k optimization parameters: begin: " << koptBegin
          << ", kopt end: " << koptEnd << ", step: " << koptStep << endl;

  return true;
}

unsigned int ReliefSeqController::GetKmax() {
  map<ClassLevel, vector<unsigned int> > classIdxMap =
          dataset->GetClassIndexes();
  map<ClassLevel, vector<unsigned int> >::const_iterator mit =
          classIdxMap.begin();
  unsigned int minClassCount = mit->second.size();
  ++mit;
  for(; mit != classIdxMap.end(); ++mit) {
    if(mit->second.size() < minClassCount) {
      minClassCount = mit->second.size();
    }
  }
  return minClassCount - 1;
}

bool ReliefSeqController::ComputeScoresKopt() {
  // set the optimization parameters from the command line parameters
  if(!SetKoptParameters()) {
    return false;
  }

  // iterate over all k's
  //vector<map<string, double> > allScores;
  vector<unsigned int> koptValues;
	bool hasNames = false;
	vector<vector<double> > allScores;
	vector<string> scoreNames;
  for(unsigned int thisK = koptBegin; thisK <= koptEnd; thisK += koptStep) {
    // run ReliefF on this k
    cout << Timestamp() << "--------------------------" << endl;
    cout << Timestamp() << "Running ReliefSeq for k=" << thisK << endl;
    reliefseqAlgorithm->SetK(thisK);
    koptValues.push_back(thisK);
    scores.clear();
		dataset->ResetNearestNeighbors();
    scores = reliefseqAlgorithm->ComputeScores();
	  sort(scores.begin(), scores.end(), scoresSortAscByName);
		vector<double> thisScores;
		AttributeScoresCIt scoresIt = scores.begin();
		for(; scoresIt != scores.end(); ++scoresIt) {
			if(!hasNames) {
				scoreNames.push_back(scoresIt->second);
			}
			thisScores.push_back(scoresIt->first);
		}
		allScores.push_back(thisScores);
	
		// I/O
    if(paramsMap.count("write-each-k-scores")) {
      stringstream filePrefix;
      filePrefix << outFilesPrefix << "." << thisK;
      WriteAttributeScores(filePrefix.str());
    }
    // PrintScores();
    hasNames = true;
  }

  // print allScores
//	for(unsigned int i=0; i < koptValues.size(); ++i) {
//		for(unsigned int j=0; j < scoreNames.size(); ++j) {
//      cout << allScores[i][j] << " ";
//    }
//    cout << endl;
//  }
  
  // pick best scores and k's for each attribute
  scores.clear();
	for(unsigned int i=0; i < scoreNames.size(); ++i) {
    string thisVar = scoreNames[i];
		unsigned int bestK = koptValues[0];
    // cout << thisVar;
		double bestScore = -1.0;
		for(unsigned int j=0; j < koptValues.size(); ++j) {
			unsigned int thisK = koptValues[j];
			double thisScore = allScores[j][i];
      if(thisScore > bestScore) {
        bestScore = thisScore;
        bestK = thisK;
      }
    }
    // cout << "\t" << bestScore << " (" << bestK << ")" << endl;
    scores.push_back(make_pair(bestScore, thisVar));
    bestKs[thisVar] = bestK;
  }

  sort(scores.begin(), scores.end(), scoresSortDesc);
  
  if(paramsMap.count("write-best-k")) {
    WriteBestKs(outFilesPrefix);
  }
  
  return true;
}

AttributeScores& ReliefSeqController::GetScores() {
  return scores;
}

string ReliefSeqController::GetAlgorithmMode() {
  return algorithmMode;
}

void ReliefSeqController::PrintAttributeScores(ofstream& outStream) {
  for(AttributeScoresCIt scoresIt = scores.begin(); scoresIt != scores.end();
          ++scoresIt) {
    outStream << (*scoresIt).first << "\t" << (*scoresIt).second << endl;
  }
}

void ReliefSeqController::PrintScores() {
  for(AttributeScoresCIt scoresIt = scores.begin();
      scoresIt != scores.end(); ++scoresIt) {
    cout << scoresIt->first << "\t" << scoresIt->second << endl;
  }
}

void ReliefSeqController::PrintBestKs() {
  for(map<string, unsigned int>::const_iterator kIt = bestKs.begin();
      kIt != bestKs.end(); ++kIt) {
    cout << kIt->first << "\t" << kIt->second << endl;
  }
}

void ReliefSeqController::WriteAttributeScores(string baseFilename) {
  string resultsFilename = baseFilename;
  ofstream outFile;
  resultsFilename = baseFilename + ".reliefseq";
  outFile.open(resultsFilename.c_str());
  if(outFile.bad()) {
    cerr << "ERROR: Could not open scores file " << resultsFilename
            << "for writing" << endl;
    exit(1);
  }
  cout << Timestamp()
          << "Writing reliefseq scores to [" + resultsFilename + "]" << endl;
  PrintAttributeScores(outFile);
  outFile.close();
}

void ReliefSeqController::WriteBestKs(string baseFilename) {
  string resultsFilename = baseFilename;
  ofstream outFile;
  resultsFilename = baseFilename + ".bestk";
  outFile.open(resultsFilename.c_str());
  if(outFile.bad()) {
    cerr << "ERROR: Could not open scores file " << resultsFilename
            << "for writing" << endl;
    exit(1);
  }
  cout << Timestamp()
          << "Writing reliefseq best k's to [" + resultsFilename + "]" << endl;
  for(map<string, unsigned int>::const_iterator kIt = bestKs.begin();
    kIt != bestKs.end(); ++kIt) {
    outFile << kIt->first << "\t" << kIt->second << endl;
  }

  outFile.close();
}

bool ReliefSeqController::RunReliefF() {

  scores = reliefseqAlgorithm->ComputeScores();
  if(scores.size() == 0) {
    cerr << "ERROR: RunReliefF: No scores computed" << endl;
    return false;
  }

  if(algorithmMode == "reliefseq") {
    cout << Timestamp() << "rnaSeq skipping normalization." << endl;
    return true;
  }

  cout << Timestamp() << "Normalizing ReliefF scores to 0-1" << endl;
  pair<double, string> firstScore = scores[0];
  double minRFScore = firstScore.first;
  double maxRFScore = firstScore.first;
  AttributeScoresCIt scoresIt = scores.begin();
  for(; scoresIt != scores.end(); ++scoresIt) {
    pair<double, string> thisScore = *scoresIt;
    if(thisScore.first < minRFScore) {
      minRFScore = thisScore.first;
    }
    if(thisScore.first > maxRFScore) {
      maxRFScore = thisScore.first;
    }
  }

  if(!reliefseqAlgorithm->DoNormalize()) {
    return true;
  }

  // normalize attribute scores if necessary
  if(minRFScore == maxRFScore) {
    cout << Timestamp() << "WARNING: Relief-F min and max scores are the same. "
            << "No normalization necessary" << endl;
    return true;
  }

  AttributeScores newRFScores;
  double rfRange = maxRFScore - minRFScore;
  for(AttributeScoresIt it = scores.begin(); it != scores.end(); ++it) {
    pair<double, string> thisScore = *it;
    double key = thisScore.first;
    string val = thisScore.second;
    newRFScores.push_back(make_pair((key - minRFScore) / rfRange, val));
  }

  scores.clear();
  scores = newRFScores;

  return true;
}

bool ReliefSeqController::RemoveWorstAttributes(unsigned int numToRemove) {
  unsigned int numToRemoveAdj = numToRemove;
  unsigned int numAttr = dataset->NumAttributes();
  if((numAttr - numToRemove) < numTargetAttributes) {
    cout << Timestamp() << "WARNING: attempt to remove " << numToRemove
            << " attributes which will remove more than target "
            << "number of attributes " << numTargetAttributes << ". Adjusting"
            << endl;
    numToRemoveAdj = numAttr - numTargetAttributes;
  }
  cout << Timestamp() << "Removing " << numToRemoveAdj << " attributes" << endl;
  sort(scores.begin(), scores.end(), scoresSortAsc);
  for(unsigned int i = 0; i < numToRemoveAdj; ++i) {

    // worst score and attribute name
    pair<double, string> worst = scores[i];
    //    cout << "\t\t\t\tRemoving: "
    //            << worst.second << " (" << worst.first << ")" << endl;

    // save worst
    removedAttributes.push_back(worst);
    // remove the attribute from those under consideration
    if(!dataset->MaskRemoveVariable(worst.second)) {
      cerr << "ERROR: Could not remove worst attribute: " << worst.second
              << endl;
      return false;
    }
  }

  return true;
}

