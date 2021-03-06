/*
 * ReliefSeqCLI.cpp - Bill White - 12/20/12
 *
 * Command line interface for running ReliefF on biological SEQuence data.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <map>
#include <vector>
#include <ctime>

#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/progress.hpp>

#include "Insilico.h"
#include "Dataset.h"
#include "DgeData.h"
#include "BirdseedData.h"
#include "ReliefSeqController.h"

using namespace std;
using namespace boost;
namespace po = boost::program_options;

int main(int argc, char** argv) {

	/// command line processing variables: defaults and storage for boost
	// file names
	string configFilename = "";
	string snpsFilename = "";
	string snpsFileType = "";
	string snpExclusionFile = "";
	string numericsFilename = "";
	string dgeCountsFilename = "";
	string dgeNormsFilename = "";
	string birdseedFilename = "";
	string birdseedPhenosFilename = "";
	string birdseedSubjectsFilename = "";
	string birdseedIncludeSnpsFilename = "";
	string birdseedExcludeSnpsFilename = "";
	string altPhenotypeFilename = "";
	string outputDatasetFilename = "";
	string outputFilesPrefix = "reliefseq_default";
	string distanceMatrixFilename = "";
	string gainMatrixFilename = "";
	string titvFilename = "";
	string diagnosticLogFilename = "";
	string diagnosticLevelsCountsFilename = "";
	// ReliefF
	unsigned int k = 10;
  unsigned int koptBegin = 1;
  unsigned int koptEnd = 1;
  unsigned int koptStep = 1;
	unsigned int m = 0;
	string snpMetric = "gm";
	string snpMetricNN = "gm";
	string snpMetricWeights = "gm";
	string numMetric = "manhattan";
	string weightByDistanceMethod = "equal";
	double weightByDistanceSigma = 2.0;
	string reliefMode = "relieff";
	string reliefSeqAlgorithmMode = "snr";
	double reliefSeqAlgorithmS0 = 0.05;
	string reliefSeqSnrMode = "snr";
	string reliefSeqTstatMode = "pval";
	unsigned int reliefNumTarget = 0;
	unsigned int reliefIterNumToRemove = 0;
	unsigned int reliefIterPercentToRemove = 0;
	// numeric data parameters
	string numericTransform = "";
	unsigned int normalizeScores = 0;
	
  /// declare the supported options
  po::options_description desc("Allowed options");
  desc.add_options()
		("help", "produce help message")
		("verbose", "verbose output")
		("convert", "convert data set to data set - does not run reliefseq")
		("write-best-k", "optimize k, write best k's")
		("write-each-k-scores", "optimize k, write best scores for each k")
		(
		"config-file,c",
		po::value<string>(&configFilename),
		"read configuration options from file - command line overrides these"
		)
		(
		"snp-data,s",
		po::value<string>(&snpsFilename),
		"read SNP attributes from genotype filename: txt, ARFF, plink (map/ped, binary, raw)"
		)
		(
		"snp-file-type",
		po::value<string>(&snpsFileType),
		"Ignore file extension and use type: textwhitesp, wekaarff, plinkped, "
		"plinkbed, plinkraw, dge, birdseed"
		)
		(
		"numeric-data,n",
		po::value<string>(&numericsFilename),
		"read continuous attributes from PLINK-style covar file"
		)
		(
		"numeric-transform,X",
		po::value<string>(&numericTransform),
		"perform numeric transformation: normalize, standardize, zscore, log, sqrt, anscombe"
		)
		(
		"alternate-pheno-file,a",
		po::value<string>(&altPhenotypeFilename),
		"specifies an alternative phenotype/class label file; one value per line"
		)
		(
		"algorithm-mode,g",
		po::value<string>(&reliefMode)->default_value(reliefMode),
		"Relief algorithm mode (relieff|reliefseq)"
		)
		(
		"seq-algorithm-mode",
		po::value<string>(&reliefSeqAlgorithmMode)->default_value(reliefSeqAlgorithmMode),
		"Relief algorithm mode (snr|tstat)"
		)
		(
		"seq-snr-mode",
		po::value<string>(&reliefSeqSnrMode)->default_value(reliefSeqSnrMode),
		"Seq interaction algorithm SNR mode (snr|relieff)"
		)
		(
		"seq-tstat-mode",
		po::value<string>(&reliefSeqTstatMode)->default_value(reliefSeqTstatMode),
		"Seq interaction algorithm t-statistic mode (pval|abst|rawt)"
		)
		(
		"seq-algorithm-s0",
		po::value<double>(&reliefSeqAlgorithmS0)->default_value(reliefSeqAlgorithmS0),
		"Seq interaction algorithm s0 (0.0 <= s0 <= 1.0)"
		)
		(
		"num-target,t",
		po::value<unsigned int>(&reliefNumTarget),
		"target number of attributes to keep after backwards selection"
		)
		(
		"iter-remove-n,r",
		po::value<unsigned int>(&reliefIterNumToRemove),
		"number of attributes to remove per iteration of backwards selection"
		)
		(
		"iter-remove-percent,p",
		po::value<unsigned int>(&reliefIterPercentToRemove),
		"percentage of attributes to remove per iteration of backwards selection"
		)
		(
		"normalize-scores",
		po::value<unsigned int>(&normalizeScores)->default_value(normalizeScores),
		"normalize ReliefF scores? (0|1)"
		)
		(
		"out-dataset-filename,O",
		po::value<string > (&outputDatasetFilename),
		"write a new tab-delimited data set with Reliefseq filtered attributes"
		)
		(
		"out-files-prefix,o",
		po::value<string > (&outputFilesPrefix)->default_value(outputFilesPrefix),
		"use prefix for all output files"
		)
		(
		"snp-metric-nn,B",
		po::value<string > (&snpMetricNN)->default_value(snpMetricNN),
		"metric for determining the difference between subjects (gm|am|nca|nca6|km)"
		)
		(
		"snp-metric-weights,W",
		po::value<string > (&snpMetricWeights)->default_value(snpMetricWeights),
		"metric for determining the diff(erence) between SNPs (gm|am|nca|nca6)"
		)
		(
		"numeric-metric,N",
		po::value<string > (&numMetric)->default_value(numMetric),
		"metric for determining the difference between numeric attributes (manhattan|euclidean)"
		)
		(
		"snp-exclusion-file,x",
		po::value<string > (&snpExclusionFile),
		"file of SNP names to be excluded"
		)
		(
		"k-nearest-neighbors,k",
		po::value<unsigned int>(&k)->default_value(k),
		"set k nearest neighbors (0=optimize k)"
		)
		(
		"kopt-begin",
		po::value<unsigned int>(&koptBegin)->default_value(koptBegin),
		"optimize k starting with kopt-begin"
		)
		(
		"kopt-end",
		po::value<unsigned int>(&koptEnd)->default_value(koptEnd),
		"optimize k ending with kopt-end"
		)
		(
		"kopt-step",
		po::value<unsigned int>(&koptStep)->default_value(koptStep),
		"optimize k incrementing with kopt-step"
		)
		(
		"number-random-samples,m",
		po::value<unsigned int>(&m)->default_value(m),
		"number of random samples (0=all|1 <= n <= number of samples)"
		)
		(
		"weight-by-distance-method,b",
		po::value<string > (&weightByDistanceMethod)->default_value(weightByDistanceMethod),
		"weight-by-distance method (equal|one_over_k|exponential)"
		)
		(
		"weight-by-distance-sigma",
		po::value<double>(&weightByDistanceSigma)->default_value(weightByDistanceSigma),
		"weight by distance sigma"
		)
		(
		"diagnostic-tests,d",
		po::value<string > (&diagnosticLogFilename),
		"performs diagnostic tests and sends output to filename without running reliefseq"
		)
		(
		"diagnostic-levels-file,D",
		po::value<string > (&diagnosticLevelsCountsFilename),
		"write diagnostic attribute level counts to filename"
		)
		(
		"dge-counts-data",
		po::value<string > (&dgeCountsFilename),
		"read digital gene expression counts from text file"
		)
		(
		"dge-norm-factors",
		po::value<string > (&dgeNormsFilename),
		"read digital gene expression normalization factors from text file"
		)
		(
		"birdseed-snps-data",
		po::value<string > (&birdseedFilename),
		"read SNP data from a birdseed formatted file"
		)
		(
		"birdseed-phenos-data",
		po::value<string > (&birdseedPhenosFilename),
		"read birdseed subjects phenotypes from a text file"
		)
		(
		"birdseed-subjects-labels",
		po::value<string > (&birdseedSubjectsFilename),
		"read subject labels from filename to override names from data file"
		)
		(
		"birdseed-include-snps",
		po::value<string > (&birdseedIncludeSnpsFilename),
		"include the SNP IDs listed in the text file"
		)
		(
		"birdseed-exclude-snps",
		po::value<string > (&birdseedExcludeSnpsFilename),
		"exclude the SNP IDs listed the text file"
		)
		(
		"distance-matrix",
		po::value<string > (&distanceMatrixFilename),
		"create a distance matrix for the loaded samples and exit"
		)
		(
		"gain-matrix",
		po::value<string > (&gainMatrixFilename),
		"create a GAIN matrix for the loaded samples and exit"
		)
		(
		"dump-titv-file",
		po::value<string > (&titvFilename),
		"file for dumping SNP transition/transversion information"
		)
		;

	/// parse the command line and/or config file into a Boost variables map
	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	/// user needs help
	if(vm.count("help") || (argc == 1)) {
		cerr << desc << endl;
		exit(COMMAND_LINE_ERROR);
	}

	/// begin main program
	// ---------------------------------------------------------------------------
	cout << Timestamp() << argv[0] << " starting" << endl;
	boost::progress_timer t;

	// ---------------------------------------------------------------------------
	cout << Timestamp() << "Processing command line arguments" << endl;

	/// read config file if specified
	if(vm.count("config-file")) {
		ifstream configStream(configFilename.c_str());
		if (!configStream.is_open()) {
			cerr << "ERROR: Could not open configuration file: "
					<< configFilename << endl;
			exit(EXIT_FAILURE);
		}
		cout << Timestamp() << "Reading configuration options from: "
				<< configFilename << endl;
		po::store(po::parse_config_file(configStream, desc), vm);
		po::notify(vm);
		configStream.close();
	}

	/// determine the output data set type
	OutputDatasetType outputDatasetType = NO_OUTPUT_DATASET;
	if(outputDatasetFilename != "") {
		// determine the data set type
		string outFileExtension = GetFileExtension(outputDatasetFilename);
		if((outFileExtension == "txt") || (outFileExtension == "tab")) {
			outputDatasetType = TAB_DELIMITED_DATASET;
		} else {
			if(outFileExtension == "csv") {
				outputDatasetType = CSV_DELIMITED_DATASET;
			} else {
				if(outFileExtension == "arff") {
					outputDatasetType = ARFF_DATASET;
				} else {
					if(outFileExtension == "ped") {
						outputDatasetType = PLINK_PED_DATASET;
					}
					else {
						if(outFileExtension == "cov") {
							outputDatasetType = PLINK_COVAR_DATASET;
						}
						else {
							cerr << "Unrecognized output data set filename extension: ["
										<< outFileExtension << "]" << endl;
							exit(COMMAND_LINE_ERROR);
						}
					}
				}
			}
		}
		cout << Timestamp() << "Writing ReliefF filtered data set to ["
						<< outputDatasetFilename
						<< "]" << endl;
	}

  /// determine the analysis type
  cout << Timestamp() << "Determining analysis type" << endl;
  AnalysisType analysisType = NO_ANALYSIS;
  bool noAnalysisFound = true;
  if(noAnalysisFound &&
  		(vm.count("diagnostic-tests") || vm.count("diagnostic-levels-file"))) {
    cout << Timestamp() << "Diagnostic test requested" << endl;
    analysisType = DIAGNOSTIC_ANALYSIS;
    noAnalysisFound = false;
  }
	if(noAnalysisFound && vm.count("convert")) {
		cout << Timestamp() << "Data set conversion requested" << endl;
		analysisType = DATASET_CONVERSION;
    noAnalysisFound = false;
	}
	if(noAnalysisFound && (vm.count("snp-data") && vm.count("numeric-data"))) {
		cout << Timestamp() << "Integrated analysis requested" << endl;
		analysisType = INTEGRATED_ANALYSIS;
    noAnalysisFound = false;
	}
	if(noAnalysisFound && (vm.count("snp-data") && !vm.count("numeric-data"))) {
		cout << Timestamp() << "SNP-only analysis requested" << endl;
		analysisType = SNP_ONLY_ANALYSIS;
    noAnalysisFound = false;
	}
	if(noAnalysisFound && (!vm.count("snp-data") && (vm.count("numeric-data")))) {
		cout << Timestamp() << "Numeric-only analysis requested" << endl;
		analysisType = NUMERIC_ONLY_ANALYSIS;
    noAnalysisFound = false;
	}
	if(noAnalysisFound && vm.count("dge-counts-data")) {
		cout << Timestamp() << "DGE analysis requested" << endl;
		analysisType = DGE_ANALYSIS;
		noAnalysisFound = false;
	}
	if(noAnalysisFound && vm.count("birdseed-snps-data")) {
		cout << Timestamp() << "Birdseed SNPs analysis requested" << endl;
		analysisType = BIRDSEED_ANALYSIS;
    noAnalysisFound = false;
	}
	if(noAnalysisFound) {
		cerr << "ERROR: Could not determine the analysis to perform based on "
						<< "command line options: " << endl << desc << endl;
		exit(COMMAND_LINE_ERROR);
	}

	/// check for numerics and alternate phenotype files; if present need to
	/// match/intersect the IDs used and only load those IDs from the data set
	// -------------------------------------------------------------------------
	// added bcw 7/15/11 - number of numerics or phenotypes might not be the
	// same as the data set - only load those in the numerics/phenotype file
	// if covariate or alternate phenotype file is present, read it first to
	// get the keys for reading the data set instances
	cout << Timestamp()
					<< "Checking for numeric data and/or alternate phenotype files" << endl;
	vector<string> numericsIds;
	vector<string> phenoIds;
	if(analysisType == SNP_ONLY_ANALYSIS ||
		 analysisType == NUMERIC_ONLY_ANALYSIS ||
		 analysisType == INTEGRATED_ANALYSIS) {
		if(numericsFilename != "") {
			cout << Timestamp() << "Loading individual IDs from numeric data file: "
							<< numericsFilename << endl;
			if(!LoadNumericIds(numericsFilename, numericsIds)) {
				exit(COMMAND_LINE_ERROR);
			}
			// copy (numericsIds.begin(), numericsIds.end(), ostream_iterator<string> (cout, "\n"));
		}
		if(altPhenotypeFilename != "") {
			cout << Timestamp() << "Loading individual IDs from alternate phenotype file: "
							<< altPhenotypeFilename << endl;
			if(!LoadPhenoIds(altPhenotypeFilename, phenoIds)) {
				exit(COMMAND_LINE_ERROR);
			}
			// copy(phenoIds.begin(), phenoIds.end(), ostream_iterator<string> (cout, "\n"));
		}
	} else {
		cout << Timestamp() << "Numeric data and alternate phenotype files are not "
				"specified for this analysis type" << endl;
	}

	// -------------------------------------------------------------------------
	/// find IDs for loading from the dataset
	cout << Timestamp() << "Determining the IDs to be read from the dataset" << endl;
	vector<string> indIds;
	if(!GetMatchingIds(numericsFilename, altPhenotypeFilename,
										 numericsIds, phenoIds, indIds)) {
		cerr << "ERROR: could not get matching IDs from numeric " <<
						" and/or phenotype files." << endl;
		exit(COMMAND_LINE_ERROR);
	}
	cout << Timestamp() << indIds.size()
					<< " individual IDs read from numeric and/or phenotype file(s)"
					<< endl;

	// -------------------------------------------------------------------------
	/// load and prepare data for running reliefseq
	cout << Timestamp() << "Loading and preparing data for Reliefseq analysis" << endl;
	Dataset* ds = 0;
	DgeData* dge = 0;
	BirdseedData* birdseed = 0;
	bool datasetLoaded = false;
	switch(analysisType) {
		case SNP_ONLY_ANALYSIS:
			cout << Timestamp() << "Reading SNPs data set" << endl;
			ds = ChooseSnpsDatasetByType(snpsFilename, snpsFileType);
			datasetLoaded = ds->LoadDataset(snpsFilename, "",
																			altPhenotypeFilename, indIds);
			break;
		case NUMERIC_ONLY_ANALYSIS:
			cout << Timestamp() << "Reading numeric data set" << endl;
			ds = new Dataset();
			datasetLoaded = ds->LoadDataset("", numericsFilename,
																			altPhenotypeFilename, indIds);
			break;
		case INTEGRATED_ANALYSIS:
			cout << Timestamp() << "Reading datasets for integrated analysis" << endl;
			ds = ChooseSnpsDatasetByType(snpsFilename, snpsFileType);
			datasetLoaded = ds->LoadDataset(snpsFilename, numericsFilename,
																			altPhenotypeFilename, indIds);
			break;
		case DGE_ANALYSIS:
			cout << Timestamp() << "Reading numerics data set from digital gene "
			<< "expression (DGE) data" << endl;
			dge = new DgeData();
			if(dge->LoadData(dgeCountsFilename, dgeNormsFilename)) {
				ds = new Dataset();
				datasetLoaded = ds->LoadDataset(dge);
			}
			else {
				cerr << "ERROR: Failure to load DGE data set"
						<< endl << endl;
			}
			break;
		case BIRDSEED_ANALYSIS:
			cout << Timestamp() << "Reading SNPs data set from Birdseed-called "
			<< "data" << endl;
			birdseed = new BirdseedData();
			if(birdseed->LoadData(birdseedFilename, birdseedPhenosFilename,
					birdseedSubjectsFilename, birdseedIncludeSnpsFilename,
					birdseedExcludeSnpsFilename)) {
				ds = new Dataset();
				datasetLoaded = ds->LoadDataset(birdseed);
			}
			else {
				cerr << "ERROR: Failure to load Birdseed data set"
						<< endl << endl;
			}
			break;
		case DIAGNOSTIC_ANALYSIS:
			cout << Timestamp() << "Performing SNP diagnostics" << endl;
			if(snpsFilename == "" && birdseedFilename == "") {
				cerr << "Cannot run diagnostics without a SNP file specified with "
								<< "--snp-data or --birdseed-snps-data" << endl;
				exit(COMMAND_LINE_ERROR);
			}
			if(snpsFilename != "") {
				ds = ChooseSnpsDatasetByType(snpsFilename, snpsFileType);
				datasetLoaded = ds->LoadDataset(snpsFilename, numericsFilename,
												altPhenotypeFilename, indIds);
			}
			else {
				birdseed = new BirdseedData();
				if(birdseed->LoadData(birdseedFilename, birdseedPhenosFilename,
						birdseedSubjectsFilename, birdseedIncludeSnpsFilename,
						birdseedExcludeSnpsFilename)) {
					ds = new Dataset();
					datasetLoaded = ds->LoadDataset(birdseed);
				}
			}
			if(!datasetLoaded) {
				cerr << "ERROR: Failure to load data set for diagnostic analysis"
						<< endl << endl;
			}
			ds->RunSnpDiagnosticTests(diagnosticLogFilename);
			if(diagnosticLevelsCountsFilename != "") {
				ds->WriteLevelCounts(diagnosticLevelsCountsFilename);
			}
			cout << Timestamp() << argv[0] << " done" << endl;
			delete ds;
			return 0;
			break;
		case DATASET_CONVERSION:
			if((snpsFilename != "" || birdseedFilename != "") &&
					numericsFilename == ""
					) {
				if(snpsFilename != "") {
					ds = ChooseSnpsDatasetByType(snpsFilename,snpsFileType);
					datasetLoaded = ds->LoadDataset(snpsFilename, numericsFilename,
													altPhenotypeFilename, indIds);
				}
				else {
					birdseed = new BirdseedData();
					if(birdseed->LoadData(birdseedFilename, birdseedPhenosFilename,
							birdseedSubjectsFilename, birdseedIncludeSnpsFilename,
							birdseedExcludeSnpsFilename)) {
						ds = new Dataset();
						datasetLoaded = ds->LoadDataset(birdseed);
					}
				}
				if(!datasetLoaded) {
					cerr << "ERROR: Failure to load data set for conversion"
							<< endl << endl;
				}
				switch(outputDatasetType) {
					case TAB_DELIMITED_DATASET:
					case CSV_DELIMITED_DATASET:
					case ARFF_DATASET:
					case PLINK_PED_DATASET:
						ds->WriteNewDataset(outputDatasetFilename, outputDatasetType);
						break;
					case PLINK_BED_DATASET:
						cout << "PLINK BED output format not supported yet" << endl;
						exit(COMMAND_LINE_ERROR);
						break;
					case PLINK_COVAR_DATASET:
						ds->WriteNewDataset(outputDatasetFilename, outputDatasetType);
						break;
					case NO_OUTPUT_DATASET:
					default:
						cerr << "No output data set specified for conversion" << endl;
						exit(COMMAND_LINE_ERROR);
						break;
				}
				cout << Timestamp() << "Conversion from " << snpsFilename
						<< " to " << outputDatasetFilename << " successful" << endl;
				cout << Timestamp() << "Conversion complete" << endl;
				cout << Timestamp() << argv[0] << " done" << endl;
				delete ds;
				return 0;
			}
			else {
				cerr << "ERROR: Data set conversion with SNP data sets"
						<< endl << endl;
				exit(COMMAND_LINE_ERROR);
			}
			break;
		case NO_ANALYSIS:
			cerr << "Analysis type could not be determined" << endl;
			exit(COMMAND_LINE_ERROR);
		default:
			cerr << "INTERNAL ERROR: Undefined analysis type: "
				<< analysisType << endl;
			exit(COMMAND_LINE_ERROR);
	}

	if(!datasetLoaded) {
		cerr << "ERROR: Failure to load dataset for analysis" << endl << endl;
		exit(DATASET_LOAD_ERROR);
	}

	if(snpExclusionFile != "") {
		if(!ds->ProcessExclusionFile(snpExclusionFile)) {
			cerr << "ERROR: processing exclusion file: " << snpExclusionFile << endl;
			exit(EXIT_FAILURE);
		}
		cout << Timestamp() << ds->NumAttributes()
				<< " SNPs remain after processing exclusion file" << endl;
	}

	if(vm.count("snp-metric-nn")) {
		snpMetricNN = snpMetricNN;
	}

	if(vm.count("snp-metric-weights")) {
		snpMetric = snpMetricWeights;
	}

	if(!ds->SetDistanceMetrics(snpMetric, snpMetricNN, numMetric)) {
		cerr << "Could not set distance metrics for the data set, "
				<< "SNP nearest neighbors: " << snpMetricNN
				<< ", SNP ReliefF weight updates: " << snpMetricWeights
				<< ", Numeric: " << numMetric << endl;
		exit(COMMAND_LINE_ERROR);
	}

	/// happy lights
	switch(analysisType) {
		case SNP_ONLY_ANALYSIS:
		case DISTANCE_MATRIX_ANALYSIS:
		case DIAGNOSTIC_ANALYSIS:
			ds->PrintStats();
			break;
		case NUMERIC_ONLY_ANALYSIS:
		case INTEGRATED_ANALYSIS:
		case REGRESSION_ANALYSIS:
		case RNASEQ_ANALYSIS:
			ds->PrintNumericsStats();
			break;
		case DGE_ANALYSIS:
			dge->PrintSampleStats();
			break;
		case BIRDSEED_ANALYSIS:
			birdseed->PrintInfo();
			break;
		case DATASET_CONVERSION:
		case NO_ANALYSIS:
			break;
	}

	if(vm.count("dump-titv-file")) {
		cout << Timestamp()
				<< "Dumping SNP transition/transversion information to ["
				<< titvFilename << "]" << endl;
		ds->WriteSnpTiTvInfo(titvFilename);
	  cout << Timestamp() << "Elapsed time: " << t.elapsed() << " secs" << endl;
	  cout << Timestamp() << argv[0] << " done" << endl;
		return 0;
	}

	/// distance matrix calculation(s) do their work, then exit main()
	if(vm.count("distance-matrix") || vm.count("gain-matrix")) {
		if(vm.count("distance-matrix")) {
			double** distanceMatrix = 0;
		  /// create a distance matrix
		  vector<string> instanceIds = ds->MaskGetInstanceIds();
		  int numInstances = instanceIds.size();
		  distanceMatrix = new double*[numInstances];
		  for(int i = 0; i < numInstances; ++i) {
		    distanceMatrix[i] = new double[numInstances];
		    for(int j = 0; j < numInstances; ++j) {
		      distanceMatrix[i][j] = 0.0;
		    }
		  }
			if(ds->CalculateDistanceMatrix(distanceMatrix, distanceMatrixFilename)) {
				for(unsigned int i=0; i < ds->NumInstances(); ++i) {
					delete [] distanceMatrix[i];
				}
				delete [] distanceMatrix;
			}
			else {
				cerr << "ERROR: Could not calculate a distance matrix." << endl;
				exit(EXIT_FAILURE);
			}
		}
		if(vm.count("gain-matrix")) {
			double** gainMatrix = 0;
		  /// create an attribute interaction matrix
		  vector<unsigned int> attributeIds =
		  		ds->MaskGetAttributeIndices(DISCRETE_TYPE);
		  int numAttributes = attributeIds.size();
		  gainMatrix = new double*[numAttributes];
		  for(int i = 0; i < numAttributes; ++i) {
		    gainMatrix[i] = new double[numAttributes];
		    for(int j = 0; j < numAttributes; ++j) {
		      gainMatrix[i][j] = 0.0;
		    }
		  }
			if(ds->CalculateGainMatrix(gainMatrix, gainMatrixFilename)) {
				for(unsigned int i=0; i < ds->NumAttributes(); ++i) {
					delete [] gainMatrix[i];
				}
				delete [] gainMatrix;
			}
			else {
				cerr << "ERROR: Could not calculate a GAIN matrix." << endl;
				exit(EXIT_FAILURE);
			}
		}
	  cout << Timestamp() << "Elapsed time: " << t.elapsed() << " secs" << endl;
	  cout << Timestamp() << argv[0] << " done" << endl;
		return 0;
	}

	// do any data transformations before the analysis
	if(ds->HasNumerics() && numericTransform != "") {
		cout << Timestamp() << "Performing numeric transformation: "
				<< numericTransform << endl;
		if(numericTransform == "standardize") {
			ds->TransformNumericsStandardize();
		}
		if(numericTransform == "normalize") {
			ds->TransformNumericsNormalize();
		}
		if(numericTransform == "log") {
			ds->TransformNumericsLog();
		}
		if(numericTransform == "sqrt") {
			ds->TransformNumericsSqrt();
		}
		if(numericTransform == "anscombe") {
			ds->TransformNumericsAnscombe();
		}
		if(numericTransform == "zscore") {
			ds->TransformNumericsZScore();
		}
		//		cout << "DEBUG: writing debug_transformed.txt, transformed by: "
		//				<< numericTransform << ", and exiting."
		//				<< endl << endl;
		//		ds->WriteNewDataset("debug_transformed.txt", TAB_DELIMITED_DATASET);
		//		exit(0);
	}

	// ---------------------------------------------------------------------------
	// FINALLY! run the algorithm
	cout << Timestamp() << "Running ReliefSeq" << endl;
	ReliefSeqController rsc(ds, vm, analysisType);
	if(k == 0) {
		if(!rsc.ComputeScoresKopt()) {
			cerr << "ERROR: Failed to calculate optimum k ReliefSeq scores" << endl;
			exit(EXIT_FAILURE);
		}
	}
	else {
		if(!rsc.ComputeScores()) {
			cerr << "ERROR: Failed to calculate ReliefSeq scores" << endl;
			exit(EXIT_FAILURE);
		}
	}
	cout << Timestamp() << "ReliefSeq done" << endl;

	// ---------------------------------------------------------------------------
	// write results files
	rsc.WriteAttributeScores(outputFilesPrefix);

	/// write the ranked/filtered attributes as a new data set
	switch(outputDatasetType) {
		case TAB_DELIMITED_DATASET:
			cout << Timestamp() << "Writing new TAB file" << endl;
			ds->WriteNewDataset(outputDatasetFilename, TAB_DELIMITED_DATASET);
			break;
		case CSV_DELIMITED_DATASET:
			cout << Timestamp() << "Writing new CSV file" << endl;
			ds->WriteNewDataset(outputDatasetFilename, CSV_DELIMITED_DATASET);
			break;
		case ARFF_DATASET:
			cout << Timestamp() << "Writing new ARFF file" << endl;
			ds->WriteNewDataset(outputDatasetFilename, ARFF_DATASET);
			break;
		case PLINK_PED_DATASET:
			cout << Timestamp() << "Writing new PLINK MAP/PED files" << endl;
			ds->WriteNewDataset(outputDatasetFilename, PLINK_PED_DATASET);
			break;
		case PLINK_BED_DATASET:
		case NO_OUTPUT_DATASET:
		default:
			break;
	}

	// ---------------------------------------------------------------------------
	cout << Timestamp() << "Clean up and shutdown" << endl;
	if(ds) {
		delete ds;
	}

	cout << Timestamp() << "ReliefSeq elapsed time " << t.elapsed()
			<< " secs" << endl;
	cout << Timestamp() << argv[0] << " done" << endl;

	return 0;
}
