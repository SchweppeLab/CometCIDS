/*
   Copyright 2013 University of Washington

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include "Common.h"
#include "CometMassSpecUtils.h"
#include "CometSearch.h"
#include "CometPostAnalysis.h"
#include "CometPreprocess.h"
#include "CometWriteOut.h"
#include "CometWriteSqt.h"
#include "CometWriteTxt.h"
#include "CometWritePepXML.h"
#include "CometWritePercolator.h"
#include "Threading.h"
#include "ThreadPool.h"
#include "CometDataInternal.h"
#include "CometSearchManager.h"
#include "CometStatus.h"
#include "CometCheckForUpdates.h"
#include <sstream>

#undef PERF_DEBUG

std::vector<Query*>           g_pvQuery;
std::vector<InputFileInfo *>  g_pvInputFiles;
std::vector<double>           g_pvDIAWindows;
StaticParams                  g_staticParams;
vector<DBIndex>               g_pvDBIndex;
MassRange                     g_massRange;
map<long long, IndexProteinStruct>    g_pvProteinNames;  // for db index
Mutex                         g_pvQueryMutex;
Mutex                         g_preprocessMemoryPoolMutex;
Mutex                         g_searchMemoryPoolMutex;
CometStatus                   g_cometStatus;

/******************************************************************************
*
* Static helper functions
*
******************************************************************************/
static void GetHostName()
{
#ifdef _WIN32
   WSADATA WSAData;
   WSAStartup(MAKEWORD(1, 0), &WSAData);

   if (gethostname(g_staticParams.szHostName, SIZE_FILE) != 0)
      strcpy(g_staticParams.szHostName, "locahost");

   WSACleanup();
#else
   if (gethostname(g_staticParams.szHostName, SIZE_FILE) != 0)
      strcpy(g_staticParams.szHostName, "locahost");
#endif

   char *pStr;
   if ((pStr = strchr(g_staticParams.szHostName, '.'))!=NULL)
      *pStr = '\0';
}

static InputType GetInputType(const char *pszFileName)
{
   int iLen = (int)strlen(pszFileName);

   if (!STRCMP_IGNORE_CASE(pszFileName + iLen - 6, ".mzXML")
         || !STRCMP_IGNORE_CASE(pszFileName + iLen - 5, ".mzML")
         || !STRCMP_IGNORE_CASE(pszFileName + iLen - 9, ".mzXML.gz")
         || !STRCMP_IGNORE_CASE(pszFileName + iLen - 8, ".mzML.gz"))

   {
      return InputType_MZXML;
   }
   else if (!STRCMP_IGNORE_CASE(pszFileName + iLen - 4, ".raw"))
   {
      return InputType_RAW;
   }
   else if (!STRCMP_IGNORE_CASE(pszFileName + iLen - 4, ".ms2")
         || !STRCMP_IGNORE_CASE(pszFileName + iLen - 5, ".cms2"))
   {
      return InputType_MS2;
   }
   else if (!STRCMP_IGNORE_CASE(pszFileName + iLen - 4, ".mgf"))
   {
      return InputType_MGF;
   }

   return InputType_UNKNOWN;
}

static bool UpdateInputFile(InputFileInfo *pFileInfo)
{
   bool bUpdateBaseName = false;
   char szTmpBaseName[SIZE_FILE];

   // Make sure not set on command line OR more than 1 input file
   // Need to do this check here before g_staticParams.inputFile is set to *pFileInfo
   if (g_staticParams.inputFile.szBaseName[0] =='\0' || g_pvInputFiles.size()>1)
      bUpdateBaseName = true;
   else
      strcpy(szTmpBaseName, g_staticParams.inputFile.szBaseName);

   g_staticParams.inputFile = *pFileInfo;

   g_staticParams.inputFile.iInputType = GetInputType(g_staticParams.inputFile.szFileName);

   if (InputType_UNKNOWN == g_staticParams.inputFile.iInputType)
   {
       return false;
   }
   int iLen = (int)strlen(g_staticParams.inputFile.szFileName);

   // per request, perform quick check to validate file still exists
   // to avoid creating stub output files in these cases.
   FILE *fp;
   if ( (fp=fopen(g_staticParams.inputFile.szFileName, "r"))==NULL)
   {
      char szErrorMsg[1024];
      sprintf(szErrorMsg,  " Error - cannot read input file \"%s\".\n",
            g_staticParams.inputFile.szFileName);
      string strErrorMsg(szErrorMsg);
      g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      return false;
   }
   else
   {
      fclose(fp);
   }

   if (bUpdateBaseName) // set individual basename from input file
   {
      char *pStr;

      strcpy(g_staticParams.inputFile.szBaseName, g_staticParams.inputFile.szFileName);

      if ( (pStr = strrchr(g_staticParams.inputFile.szBaseName, '.')))
         *pStr = '\0';

      if (!STRCMP_IGNORE_CASE(g_staticParams.inputFile.szFileName + iLen - 9, ".mzXML.gz")
            || !STRCMP_IGNORE_CASE(g_staticParams.inputFile.szFileName + iLen - 8, ".mzML.gz"))
      {
         if ( (pStr = strrchr(g_staticParams.inputFile.szBaseName, '.')))
            *pStr = '\0';
      }
   }
   else
   {
      strcpy(g_staticParams.inputFile.szBaseName, szTmpBaseName);  // set basename from command line
   }

   // Create .out directory.
   if (g_staticParams.options.bOutputOutFiles)
   {
#ifdef _WIN32
      if (_mkdir(g_staticParams.inputFile.szBaseName) == -1)
      {
         errno_t err;
         _get_errno(&err);

         if (err != EEXIST)
         {
            char szErrorMsg[1024];
            sprintf(szErrorMsg,  " Error - could not create directory \"%s\".\n",
                  g_staticParams.inputFile.szBaseName);
            string strErrorMsg(szErrorMsg);
            g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
            logerr(szErrorMsg);
            return false;
         }
      }
      if (g_staticParams.options.iDecoySearch == 2)
      {
         char szDecoyDir[SIZE_FILE];
         sprintf(szDecoyDir, "%s_decoy", g_staticParams.inputFile.szBaseName);

         if (_mkdir(szDecoyDir) == -1)
         {
            errno_t err;
            _get_errno(&err);

            if (err != EEXIST)
            {
               char szErrorMsg[1024];
               sprintf(szErrorMsg,  " Error - could not create directory \"%s\".\n",  szDecoyDir);
               string strErrorMsg(szErrorMsg);
               g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
               logerr(szErrorMsg);
               return false;
            }
         }
      }
#else
      if ((mkdir(g_staticParams.inputFile.szBaseName, 0775) == -1) && (errno != EEXIST))
      {
         char szErrorMsg[1024];
         sprintf(szErrorMsg,  " Error - could not create directory \"%s\".\n", g_staticParams.inputFile.szBaseName);
         string strErrorMsg(szErrorMsg);
         g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
         logerr(szErrorMsg);
         return false;
      }
      if (g_staticParams.options.iDecoySearch == 2)
      {
         char szDecoyDir[1024];
         sprintf(szDecoyDir, "%s_decoy", g_staticParams.inputFile.szBaseName);

         if ((mkdir(szDecoyDir , 0775) == -1) && (errno != EEXIST))
         {
            char szErrorMsg[1024];
            sprintf(szErrorMsg,  " Error - could not create directory \"%s\".\n",  szDecoyDir);
            string strErrorMsg(szErrorMsg);
            g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
            logerr(szErrorMsg);
            return false;
         }
      }
#endif
   }

   return true;
}

static void SetMSLevelFilter(MSReader &mstReader)
{
   vector<MSSpectrumType> msLevel;

   if (g_staticParams.options.iMSLevel == 3)
      msLevel.push_back(MS3);
   else
      msLevel.push_back(MS2);

   mstReader.setFilter(msLevel);
}

// Allocate memory for the _pResults struct for each g_pvQuery entry.
static bool AllocateResultsMem()
{
   for (std::vector<Query*>::iterator it = g_pvQuery.begin(); it != g_pvQuery.end(); ++it)
   {
      Query* pQuery = *it;

      try
      {
         pQuery->_pResults = new Results[g_staticParams.options.iNumStored];
      }
      catch (std::bad_alloc& ba)
      {
         char szErrorMsg[256];
         sprintf(szErrorMsg, " Error - new(_pResults[]). bad_alloc: %s.\n", ba.what());
         string strErrorMsg(szErrorMsg);
         g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
         logerr(szErrorMsg);
         return false;
      }

      if (g_staticParams.options.iDecoySearch==2)
      {
         try
         {
            pQuery->_pDecoys = new Results[g_staticParams.options.iNumStored];
         }
         catch (std::bad_alloc& ba)
         {
            char szErrorMsg[256];
            sprintf(szErrorMsg, " Error - new(_pDecoys[]). bad_alloc: %s.\n", ba.what());
            string strErrorMsg(szErrorMsg);
            g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
            logerr(szErrorMsg);
            return false;
         }
      }

      for (int j=0; j<g_staticParams.options.iNumStored; j++)
      {
         pQuery->_pResults[j].dPepMass = 0.0;
         pQuery->_pResults[j].dExpect = 0.0;
         pQuery->_pResults[j].fScoreSp = 0.0;
         pQuery->_pResults[j].fXcorr = XCORR_CUTOFF;
         pQuery->_pResults[j].iLenPeptide = 0;
         pQuery->_pResults[j].iRankSp = 0;
         pQuery->_pResults[j].iMatchedIons = 0;
         pQuery->_pResults[j].iTotalIons = 0;
         pQuery->_pResults[j].szPeptide[0] = '\0';
         pQuery->_pResults[j].strSingleSearchProtein = "";
         pQuery->_pResults[j].pWhichProtein.clear();
         pQuery->_pResults[j].cPeffOrigResidue = '\0';
         pQuery->_pResults[j].iPeffOrigResiduePosition = -9;

         if (g_staticParams.options.iDecoySearch)
            pQuery->_pResults[j].pWhichDecoyProtein.clear();

         if (g_staticParams.options.iDecoySearch==2)
         {
            pQuery->_pDecoys[j].dPepMass = 0.0;
            pQuery->_pDecoys[j].dExpect = 0.0;
            pQuery->_pDecoys[j].fScoreSp = 0.0;
            pQuery->_pDecoys[j].fXcorr = XCORR_CUTOFF;
            pQuery->_pDecoys[j].iLenPeptide = 0;
            pQuery->_pDecoys[j].iRankSp = 0;
            pQuery->_pDecoys[j].iMatchedIons = 0;
            pQuery->_pDecoys[j].iTotalIons = 0;
            pQuery->_pDecoys[j].szPeptide[0] = '\0';
            pQuery->_pDecoys[j].strSingleSearchProtein = "";
            pQuery->_pDecoys[j].cPeffOrigResidue = '\0';
            pQuery->_pDecoys[j].iPeffOrigResiduePosition = -9;
         }
      }
   }

   return true;
}

static bool compareByPeptideMass(Query const* a, Query const* b)
{
   return (a->_pepMassInfo.dExpPepMass < b->_pepMassInfo.dExpPepMass);
}

static bool compareByMangoIndex(Query const* a, Query const* b)
{
   return (a->dMangoIndex < b->dMangoIndex);
}

static bool compareByScanNumber(Query const* a, Query const* b)
{
   // sort by charge state if same scan number
   if (a->_spectrumInfoInternal.iScanNumber == b->_spectrumInfoInternal.iScanNumber)
      return (a->_spectrumInfoInternal.iChargeState < b->_spectrumInfoInternal.iChargeState);
   return (a->_spectrumInfoInternal.iScanNumber < b->_spectrumInfoInternal.iScanNumber);
}

static void CalcRunTime(time_t tStartTime)
{
   char szOutFileTimeString[512];
   time_t tEndTime;
   int iTmp;

   time(&tEndTime);

   int iElapseTime=(int)difftime(tEndTime, tStartTime);

   // Print out header/search info.
   sprintf(szOutFileTimeString, "%s,", g_staticParams.szDate);
   if ( (iTmp = (int)(iElapseTime/3600) )>0)
      sprintf(szOutFileTimeString+strlen(szOutFileTimeString), " %d hr.", iTmp);
   if ( (iTmp = (int)((iElapseTime-(int)(iElapseTime/3600)*3600)/60) )>0)
      sprintf(szOutFileTimeString+strlen(szOutFileTimeString), " %d min.", iTmp);
   if ( (iTmp = (int)((iElapseTime-((int)(iElapseTime/3600))*3600)%60) )>0)
      sprintf(szOutFileTimeString+strlen(szOutFileTimeString), " %d sec.", iTmp);
   if (iElapseTime == 0)
      sprintf(szOutFileTimeString+strlen(szOutFileTimeString), " 0 sec.");
   sprintf(szOutFileTimeString+strlen(szOutFileTimeString), " on %s", g_staticParams.szHostName);

   g_staticParams.iElapseTime = iElapseTime;
   strncpy(g_staticParams.szOutFileTimeString, szOutFileTimeString, 256);
   g_staticParams.szOutFileTimeString[255]='\0';
}

// for .out files
static void PrintOutfileHeader()
{
   // print parameters

   char szIsotope[16];
   char szPeak[16];

   sprintf(g_staticParams.szIonSeries, "ion series ABCXYZ nl: %d%d%d%d%d%d %d",
         g_staticParams.ionInformation.iIonVal[ION_SERIES_A],
         g_staticParams.ionInformation.iIonVal[ION_SERIES_B],
         g_staticParams.ionInformation.iIonVal[ION_SERIES_C],
         g_staticParams.ionInformation.iIonVal[ION_SERIES_X],
         g_staticParams.ionInformation.iIonVal[ION_SERIES_Y],
         g_staticParams.ionInformation.iIonVal[ION_SERIES_Z],
         g_staticParams.ionInformation.bUseNeutralLoss);

   char szUnits[8];
   char szDecoy[20];
   char szReadingFrame[20];
   char szRemovePrecursor[20];

   if (g_staticParams.tolerances.iMassToleranceUnits==0)
      strcpy(szUnits, " AMU");
   else if (g_staticParams.tolerances.iMassToleranceUnits==1)
      strcpy(szUnits, " MMU");
   else
      strcpy(szUnits, " PPM");

   if (g_staticParams.options.iDecoySearch)
      sprintf(szDecoy, " DECOY%d", g_staticParams.options.iDecoySearch);
   else
      szDecoy[0]=0;

   if (g_staticParams.options.iRemovePrecursor)
      sprintf(szRemovePrecursor, " REMOVEPREC%d", g_staticParams.options.iRemovePrecursor);
   else
      szRemovePrecursor[0]=0;

   if (g_staticParams.options.iWhichReadingFrame)
      sprintf(szReadingFrame, " FRAME%d", g_staticParams.options.iWhichReadingFrame);
   else
      szReadingFrame[0]=0;

   if (g_staticParams.tolerances.iIsotopeError > 0)
      sprintf(szIsotope, "ISOTOPE%d", g_staticParams.tolerances.iIsotopeError);

   szPeak[0]='\0';
   if (g_staticParams.ionInformation.iTheoreticalFragmentIons==1)
      strcpy(szPeak, "PEAK1");

   sprintf(g_staticParams.szDisplayLine, "display top %d, %s%s%s%s%s%s%s%s",
         g_staticParams.options.iNumPeptideOutputLines,
         szRemovePrecursor,
         szReadingFrame,
         szPeak,
         szUnits,
         (g_staticParams.tolerances.iMassToleranceType==0?" MH+":" m/z"),
         szIsotope,
         szDecoy,
         (g_staticParams.options.bClipNtermMet?" CLIPMET":"") );
}

static bool ValidateOutputFormat()
{
   if (!g_staticParams.options.bOutputSqtStream
         && !g_staticParams.options.bOutputSqtFile
         && !g_staticParams.options.bOutputTxtFile
         && !g_staticParams.options.bOutputPepXMLFile
         && !g_staticParams.options.bOutputPercolatorFile
         && !g_staticParams.options.bOutputOutFiles)
   {
      string strError = " Please specify at least one output format.";

      g_cometStatus.SetStatus(CometResult_Failed, strError);
      string strErrorFormat = strError + "\n";
      logerr(strErrorFormat.c_str());
      return false;
   }

   return true;
}

static bool ValidateSequenceDatabaseFile()
{
   FILE *fpcheck;

   // Quick sanity check to make sure sequence db file is present before spending
   // time reading & processing spectra and then reporting this error.
   if ((fpcheck=fopen(g_staticParams.databaseInfo.szDatabase, "r")) == NULL)
   {
      char szErrorMsg[1024];
      sprintf(szErrorMsg, " Error - cannot read database file \"%s\".\n Check that the file exists and is readable.\n",
            g_staticParams.databaseInfo.szDatabase);
      string strErrorMsg(szErrorMsg);
      g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      return false;
   }

   // At this point, check extension to set whether index database or not
   if (!strcmp(g_staticParams.databaseInfo.szDatabase+strlen(g_staticParams.databaseInfo.szDatabase)-4, ".idx"))
      g_staticParams.bIndexDb = 1;

   if (g_staticParams.options.bCreateIndex && g_staticParams.bIndexDb)
   {
      char szErrorMsg[1024];
      sprintf(szErrorMsg, " Error - input database already indexed: \"%s\".\n", g_staticParams.databaseInfo.szDatabase);
      string strErrorMsg(szErrorMsg);
      g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      return false;
   }

   fclose(fpcheck);
   return true;
}

static bool ValidateScanRange()
{
   if (g_staticParams.options.scanRange.iEnd < g_staticParams.options.scanRange.iStart && g_staticParams.options.scanRange.iEnd != 0)
   {
      char szErrorMsg[256];
      sprintf(szErrorMsg, " Error - start scan is %d but end scan is %d.\n The end scan must be >= to the start scan.\n",
            g_staticParams.options.scanRange.iStart,
            g_staticParams.options.scanRange.iEnd);
      string strErrorMsg(szErrorMsg);
      g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      return false;
   }

   return true;
}

static bool ValidatePeptideLengthRange()
{
   if (g_staticParams.options.peptideLengthRange.iEnd < g_staticParams.options.peptideLengthRange.iStart && g_staticParams.options.peptideLengthRange.iEnd != 0)
   {
      char szErrorMsg[256];
      sprintf(szErrorMsg, " Error - peptide length range set as %d to %d.\n The maximum length must be >= to the minimum length.\n",
            g_staticParams.options.peptideLengthRange.iStart,
            g_staticParams.options.peptideLengthRange.iEnd);
      string strErrorMsg(szErrorMsg);
      g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      return false;
   }

   return true;
}

/******************************************************************************
*
* CometSearchManager class implementation.
*
******************************************************************************/

CometSearchManager::CometSearchManager() :
    singleSearchInitializationComplete(false),
    singleSearchThreadCount(1)
{
   // Initialize the mutexes we'll use to protect global data.
   Threading::CreateMutex(&g_pvQueryMutex);

   // Initialize the mutex we'll use to protect the preprocess memory pool
   Threading::CreateMutex(&g_preprocessMemoryPoolMutex);

   // Initialize the mutex we'll use to protect the search memory pool
   Threading::CreateMutex(&g_searchMemoryPoolMutex);

   // Initialize the Comet version
   SetParam("# comet_version ", comet_version, comet_version);
}

CometSearchManager::~CometSearchManager()
{
   // Destroy the mutex we used to protect g_pvQuery.
   Threading::DestroyMutex(g_pvQueryMutex);

   // Destroy the mutex we used to protect the preprocess memory pool
   Threading::DestroyMutex(g_preprocessMemoryPoolMutex);

   // Destroy the mutex we used to protect the search memory pool
   Threading::DestroyMutex(g_searchMemoryPoolMutex);

   //std::vector calls destructor of every element it contains when clear() is called
   g_pvInputFiles.clear();

   _mapStaticParams.clear();
}

void parseCsvString(string in, vector<int> &out) {
    out.clear();
    std::stringstream ss(in);

    for (int i; ss >> i;) {
        out.push_back(i);
        if (ss.peek() == ',') {
            ss.ignore();
        }
    }
}

bool CometSearchManager::InitializeStaticParams()
{
   int iIntData;
   double dDoubleData;
   string strData;
   IntRange intRangeData;
   DoubleRange doubleRangeData;

   if (GetParamValue("database_name", strData))
      strcpy(g_staticParams.databaseInfo.szDatabase, strData.c_str());

   if (GetParamValue("decoy_prefix", strData))
      strcpy(g_staticParams.szDecoyPrefix, strData.c_str());

   if (GetParamValue("output_suffix", strData))
      strcpy(g_staticParams.szOutputSuffix, strData.c_str());

   if (GetParamValue("dia_windows_file", strData))
      strcpy(g_staticParams.szDIAWindowsFile, strData.c_str());

   if (GetParamValue("peff_obo", strData))
      strcpy(g_staticParams.peffInfo.szPeffOBO, strData.c_str());

   GetParamValue("peff_format", g_staticParams.peffInfo.iPeffSearch);

   GetParamValue("mass_offsets", g_staticParams.vectorMassOffsets);

   GetParamValue("xcorr_processing_offset", g_staticParams.iXcorrProcessingOffset);

   GetParamValue("nucleotide_reading_frame", g_staticParams.options.iWhichReadingFrame);

   GetParamValue("mass_type_parent", g_staticParams.massUtility.bMonoMassesParent);

   GetParamValue("mass_type_fragment", g_staticParams.massUtility.bMonoMassesFragment);

   GetParamValue("show_fragment_ions", g_staticParams.options.bShowFragmentIons);

   GetParamValue("num_threads", g_staticParams.options.iNumThreads);

   GetParamValue("clip_nterm_methionine", g_staticParams.options.bClipNtermMet);

   GetParamValue("theoretical_fragment_ions", g_staticParams.ionInformation.iTheoreticalFragmentIons);
   if ((g_staticParams.ionInformation.iTheoreticalFragmentIons < 0)
         || (g_staticParams.ionInformation.iTheoreticalFragmentIons > 1))
   {
      g_staticParams.ionInformation.iTheoreticalFragmentIons = 0;
   }

   GetParamValue("use_A_ions", g_staticParams.ionInformation.iIonVal[ION_SERIES_A]);

   GetParamValue("use_B_ions", g_staticParams.ionInformation.iIonVal[ION_SERIES_B]);

   GetParamValue("use_C_ions", g_staticParams.ionInformation.iIonVal[ION_SERIES_C]);

   GetParamValue("use_X_ions", g_staticParams.ionInformation.iIonVal[ION_SERIES_X]);

   GetParamValue("use_Y_ions", g_staticParams.ionInformation.iIonVal[ION_SERIES_Y]);

   GetParamValue("use_Z_ions", g_staticParams.ionInformation.iIonVal[ION_SERIES_Z]);

   GetParamValue("use_NL_ions", g_staticParams.ionInformation.bUseNeutralLoss);

   string waterLossesOption;
   if(GetParamValue("water_losses", waterLossesOption)) {
       parseCsvString(waterLossesOption, g_staticParams.ionInformation.waterLosses);
   }
   else {
       // Default value.
       g_staticParams.ionInformation.waterLosses.clear();
       g_staticParams.ionInformation.waterLosses.push_back(1);
   }

   string ammoniaLossesOption;
   if(GetParamValue("ammonia_losses", ammoniaLossesOption)) {
       parseCsvString(ammoniaLossesOption, g_staticParams.ionInformation.ammoniaLosses);
   }
   else {
       // Default value.
       g_staticParams.ionInformation.ammoniaLosses.clear();
       g_staticParams.ionInformation.ammoniaLosses.push_back(1);
   }


   GetParamValue("variable_mod01", g_staticParams.variableModParameters.varModList[VMOD_1_INDEX]);

   GetParamValue("variable_mod02", g_staticParams.variableModParameters.varModList[VMOD_2_INDEX]);

   GetParamValue("variable_mod03", g_staticParams.variableModParameters.varModList[VMOD_3_INDEX]);

   GetParamValue("variable_mod04", g_staticParams.variableModParameters.varModList[VMOD_4_INDEX]);

   GetParamValue("variable_mod05", g_staticParams.variableModParameters.varModList[VMOD_5_INDEX]);

   GetParamValue("variable_mod06", g_staticParams.variableModParameters.varModList[VMOD_6_INDEX]);

   GetParamValue("variable_mod07", g_staticParams.variableModParameters.varModList[VMOD_7_INDEX]);

   GetParamValue("variable_mod08", g_staticParams.variableModParameters.varModList[VMOD_8_INDEX]);

   GetParamValue("variable_mod09", g_staticParams.variableModParameters.varModList[VMOD_9_INDEX]);

   if (GetParamValue("max_variable_mods_in_peptide", iIntData))
   {
      if (iIntData >= 0)
         g_staticParams.variableModParameters.iMaxVarModPerPeptide = iIntData;
   }

   GetParamValue("require_variable_mod", g_staticParams.variableModParameters.bRequireVarMod);

   GetParamValue("fragment_bin_tol", g_staticParams.tolerances.dFragmentBinSize);
   if (g_staticParams.tolerances.dFragmentBinSize < 0.01)
      g_staticParams.tolerances.dFragmentBinSize = 0.01;

   GetParamValue("fragment_bin_offset", g_staticParams.tolerances.dFragmentBinStartOffset);

   GetParamValue("peptide_mass_tolerance", g_staticParams.tolerances.dInputTolerance);

   GetParamValue("precursor_tolerance_type", g_staticParams.tolerances.iMassToleranceType);
   if ((g_staticParams.tolerances.iMassToleranceType < 0)
         || (g_staticParams.tolerances.iMassToleranceType > 1))
   {
      g_staticParams.tolerances.iMassToleranceType = 0;
   }

   GetParamValue("peptide_mass_units", g_staticParams.tolerances.iMassToleranceUnits);
   if ((g_staticParams.tolerances.iMassToleranceUnits < 0)
         || (g_staticParams.tolerances.iMassToleranceUnits > 2))
   {
      g_staticParams.tolerances.iMassToleranceUnits = 0;  // 0=amu, 1=mmu, 2=ppm
   }

   GetParamValue("isotope_error", g_staticParams.tolerances.iIsotopeError);
   if ((g_staticParams.tolerances.iIsotopeError < 0)
         || (g_staticParams.tolerances.iIsotopeError > 7))
   {
      g_staticParams.tolerances.iIsotopeError = 0;
   }

   GetParamValue("num_output_lines", g_staticParams.options.iNumPeptideOutputLines);

   GetParamValue("num_results", g_staticParams.options.iNumStored);

   GetParamValue("remove_precursor_peak", g_staticParams.options.iRemovePrecursor);

   GetParamValue("remove_precursor_tolerance", g_staticParams.options.dRemovePrecursorTol);

   if (GetParamValue("clear_mz_range", doubleRangeData))
   {
      if ((doubleRangeData.dEnd >= doubleRangeData.dStart) && (doubleRangeData.dStart >= 0.0))
      {
         g_staticParams.options.clearMzRange.dStart = doubleRangeData.dStart;
         g_staticParams.options.clearMzRange.dEnd = doubleRangeData.dEnd;
      }
   }


   GetParamValue("print_expect_score", g_staticParams.options.bPrintExpectScore);

   GetParamValue("output_sqtstream", g_staticParams.options.bOutputSqtStream);

   GetParamValue("output_sqtfile", g_staticParams.options.bOutputSqtFile);

   GetParamValue("output_txtfile", g_staticParams.options.bOutputTxtFile);

   GetParamValue("output_pepxmlfile", g_staticParams.options.bOutputPepXMLFile);

   GetParamValue("output_percolatorfile", g_staticParams.options.bOutputPercolatorFile);

   GetParamValue("output_outfiles", g_staticParams.options.bOutputOutFiles);

   GetParamValue("skip_researching", g_staticParams.options.bSkipAlreadyDone);

// GetParamValue("skip_updatecheck", g_staticParams.options.bSkipUpdateCheck);

   GetParamValue("mango_search", g_staticParams.options.bMango);

   GetParamValue("create_index", g_staticParams.options.bCreateIndex);

   GetParamValue("max_iterations", g_staticParams.options.lMaxIterations);

   GetParamValue("peff_verbose_output", g_staticParams.options.bVerboseOutput);

   GetParamValue("add_Cterm_peptide", g_staticParams.staticModifications.dAddCterminusPeptide);

   GetParamValue("add_Nterm_peptide", g_staticParams.staticModifications.dAddNterminusPeptide);

   GetParamValue("add_Cterm_protein", g_staticParams.staticModifications.dAddCterminusProtein);

   GetParamValue("add_Nterm_protein", g_staticParams.staticModifications.dAddNterminusProtein);

   if (GetParamValue("add_G_glycine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'G'] = dDoubleData;

   if (GetParamValue("add_A_alanine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'A'] = dDoubleData;

   if (GetParamValue("add_S_serine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'S'] = dDoubleData;

   if (GetParamValue("add_P_proline", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'P'] = dDoubleData;

   if (GetParamValue("add_V_valine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'V'] = dDoubleData;

   if (GetParamValue("add_T_threonine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'T'] = dDoubleData;

   if (GetParamValue("add_C_cysteine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'C'] = dDoubleData;

   if (GetParamValue("add_L_leucine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'L'] = dDoubleData;

   if (GetParamValue("add_I_isoleucine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'I'] = dDoubleData;

   if (GetParamValue("add_N_asparagine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'N'] = dDoubleData;

   if (GetParamValue("add_O_ornithine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'O'] = dDoubleData;

   if (GetParamValue("add_D_aspartic_acid", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'D'] = dDoubleData;

   if (GetParamValue("add_Q_glutamine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'Q'] = dDoubleData;

   if (GetParamValue("add_K_lysine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'K'] = dDoubleData;

   if (GetParamValue("add_E_glutamic_acid", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'E'] = dDoubleData;

   if (GetParamValue("add_M_methionine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'M'] = dDoubleData;

   if (GetParamValue("add_H_histidine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'H'] = dDoubleData;

   if (GetParamValue("add_F_phenylalanine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'F'] = dDoubleData;

   if (GetParamValue("add_R_arginine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'R'] = dDoubleData;

   if (GetParamValue("add_Y_tyrosine", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'Y'] = dDoubleData;

   if (GetParamValue("add_W_tryptophan", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'W'] = dDoubleData;

   if (GetParamValue("add_B_user_amino_acid", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'B'] = dDoubleData;

   if (GetParamValue("add_J_user_amino_acid", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'J'] = dDoubleData;

   if (GetParamValue("add_U_user_amino_acid", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'U'] = dDoubleData;

   if (GetParamValue("add_X_user_amino_acid", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'X'] = dDoubleData;

   if (GetParamValue("add_Z_user_amino_acid", dDoubleData))
      g_staticParams.staticModifications.pdStaticMods[(int)'Z'] = dDoubleData;

   GetParamValue("num_enzyme_termini", g_staticParams.options.iEnzymeTermini);
   if ((g_staticParams.options.iEnzymeTermini != 1)
         && (g_staticParams.options.iEnzymeTermini != 8)
         && (g_staticParams.options.iEnzymeTermini != 9))
   {
      g_staticParams.options.iEnzymeTermini = 2;
   }

   if (GetParamValue("scan_range", intRangeData))
   {
      if ((intRangeData.iEnd >= intRangeData.iStart) && (intRangeData.iStart > 0))
      {
         g_staticParams.options.scanRange.iStart = intRangeData.iStart;
         g_staticParams.options.scanRange.iEnd = intRangeData.iEnd;
      }
   }

   if (GetParamValue("peptide_length_range", intRangeData))
   {
      if ((intRangeData.iEnd >= intRangeData.iStart) && (intRangeData.iStart > 0))
      {
         g_staticParams.options.peptideLengthRange.iStart = intRangeData.iStart;
         g_staticParams.options.peptideLengthRange.iEnd = intRangeData.iEnd;
      }
   }

   if (GetParamValue("spectrum_batch_size", iIntData))
   {
      if (iIntData >= 0)
         g_staticParams.options.iSpectrumBatchSize = iIntData;
   }

   iIntData = 0;
   if (GetParamValue("minimum_peaks", iIntData))
   {
      if (iIntData >= 0)
         g_staticParams.options.iMinPeaks = iIntData;
   }

   if (GetParamValue("override_charge", iIntData))
   {
      if (iIntData > 0)
         g_staticParams.options.bOverrideCharge = iIntData;
   }

   if (GetParamValue("correct_mass", iIntData))
   {
      if (iIntData > 0)
         g_staticParams.options.bCorrectMass = iIntData;
   }

   if (GetParamValue("equal_I_and_L", iIntData))
   {
      g_staticParams.options.bTreatSameIL = iIntData;
   }

   if (GetParamValue("precursor_charge", intRangeData))
   {
      if ((intRangeData.iStart > 0) && (intRangeData.iEnd >= intRangeData.iStart))
      {
         g_staticParams.options.iStartCharge = intRangeData.iStart;
         g_staticParams.options.iEndCharge = intRangeData.iEnd;
      }
   }

   iIntData = 0;
   if (GetParamValue("max_fragment_charge", iIntData))
   {
      if (iIntData > MAX_FRAGMENT_CHARGE)
         iIntData = MAX_FRAGMENT_CHARGE;

      if (iIntData > 0)
         g_staticParams.options.iMaxFragmentCharge = iIntData;

      // else will go to default value (3)
   }

   iIntData = 0;
   if (GetParamValue("max_precursor_charge", iIntData))
   {
      if (iIntData > MAX_PRECURSOR_CHARGE)
         iIntData = MAX_PRECURSOR_CHARGE;

      if (iIntData > 0)
         g_staticParams.options.iMaxPrecursorCharge = iIntData;

      // else will go to default value (6)
   }

   if (GetParamValue("digest_mass_range", doubleRangeData))
   {
      if ((doubleRangeData.dEnd >= doubleRangeData.dStart) && (doubleRangeData.dStart >= 0.0))
      {
         g_staticParams.options.dPeptideMassLow = doubleRangeData.dStart;
         g_staticParams.options.dPeptideMassHigh = doubleRangeData.dEnd;
      }
   }

   if (GetParamValue("ms_level", iIntData))
   {
      if (iIntData == 3)
         g_staticParams.options.iMSLevel = 3;

      // else will go to default value (2)
   }

   if (GetParamValue("activation_method", strData))
      strcpy(g_staticParams.options.szActivationMethod, strData.c_str());

   GetParamValue("minimum_intensity", g_staticParams.options.dMinIntensity);
   if (g_staticParams.options.dMinIntensity < 0.0)
      g_staticParams.options.dMinIntensity = 0.0;

   GetParamValue("decoy_search", g_staticParams.options.iDecoySearch);
   if ((g_staticParams.options.iDecoySearch < 0) || (g_staticParams.options.iDecoySearch > 2))
      g_staticParams.options.iDecoySearch = 0;

   // Set dInverseBinWidth to its inverse in order to use a multiply instead of divide in BIN macro.
   // Safe to divide by dFragmentBinSize because of check earlier where minimum value is 0.01.
   g_staticParams.dInverseBinWidth = 1.0 /g_staticParams.tolerances.dFragmentBinSize;
   g_staticParams.dOneMinusBinOffset = 1.0 - g_staticParams.tolerances.dFragmentBinStartOffset;

   // Set masses to either average or monoisotopic.
   CometMassSpecUtils::AssignMass(g_staticParams.massUtility.pdAAMassParent,
                                  g_staticParams.massUtility.bMonoMassesParent,
                                  &g_staticParams.massUtility.dOH2parent);

   CometMassSpecUtils::AssignMass(g_staticParams.massUtility.pdAAMassFragment,
                                  g_staticParams.massUtility.bMonoMassesFragment,
                                  &g_staticParams.massUtility.dOH2fragment);

   g_staticParams.massUtility.dCO = g_staticParams.massUtility.pdAAMassFragment[(int)'c']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'o'];

   g_staticParams.massUtility.dH2O = g_staticParams.massUtility.pdAAMassFragment[(int)'h']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'h']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'o'];

   g_staticParams.massUtility.dNH3 = g_staticParams.massUtility.pdAAMassFragment[(int)'n']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'h']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'h']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'h'];

   g_staticParams.massUtility.dNH2 = g_staticParams.massUtility.pdAAMassFragment[(int)'n']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'h']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'h'];

   g_staticParams.massUtility.dCOminusH2 = g_staticParams.massUtility.dCO
            - g_staticParams.massUtility.pdAAMassFragment[(int)'h']
            - g_staticParams.massUtility.pdAAMassFragment[(int)'h'];

   GetHostName();

   // If # threads not specified, poll system to get # threads to launch.
   if (g_staticParams.options.iNumThreads <= 0)
   {
      int iNumCPUCores;
#ifdef _WIN32
      SYSTEM_INFO sysinfo;
      GetSystemInfo( &sysinfo );
      iNumCPUCores = sysinfo.dwNumberOfProcessors;

      // if user specifies a negative # threads, subtract this from # cores
      if (g_staticParams.options.iNumThreads < 0)
         g_staticParams.options.iNumThreads = iNumCPUCores + g_staticParams.options.iNumThreads;
      else
         g_staticParams.options.iNumThreads = iNumCPUCores;
#else
      iNumCPUCores= sysconf( _SC_NPROCESSORS_ONLN );

      if (g_staticParams.options.iNumThreads < 0)
         g_staticParams.options.iNumThreads = iNumCPUCores + g_staticParams.options.iNumThreads;
      else
         g_staticParams.options.iNumThreads = iNumCPUCores;

      // if set, use the environment variable NSLOTS which is defined in the qsub command
      const char * nSlots = ::getenv("NSLOTS");
      if (nSlots != NULL)
      {
         int detectedThreads = atoi(nSlots);
         if (detectedThreads > 0)
            g_staticParams.options.iNumThreads = detectedThreads;
      }
#endif
      if (g_staticParams.options.iNumThreads < 0)
      {
         g_staticParams.options.iNumThreads = 4;
         logout(" Setting number of threads to 4");
      }

      if (g_staticParams.options.iNumThreads > MAX_THREADS)
      {
         char szOut[64];
         g_staticParams.options.iNumThreads = MAX_THREADS;
         sprintf(szOut, " Setting number of threads to %d", MAX_THREADS);
         logout(szOut);
      }
   }

   // Set masses to either average or monoisotopic.
   CometMassSpecUtils::AssignMass(g_staticParams.massUtility.pdAAMassParent,
                                  g_staticParams.massUtility.bMonoMassesParent,
                                  &g_staticParams.massUtility.dOH2parent);

   CometMassSpecUtils::AssignMass(g_staticParams.massUtility.pdAAMassFragment,
                                  g_staticParams.massUtility.bMonoMassesFragment,
                                  &g_staticParams.massUtility.dOH2fragment);

   g_staticParams.massUtility.dCO = g_staticParams.massUtility.pdAAMassFragment[(int)'c']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'o'];

   g_staticParams.massUtility.dH2O = g_staticParams.massUtility.pdAAMassFragment[(int)'h']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'h']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'o'];

   g_staticParams.massUtility.dNH3 = g_staticParams.massUtility.pdAAMassFragment[(int)'n']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'h']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'h']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'h'];

   g_staticParams.massUtility.dNH2 = g_staticParams.massUtility.pdAAMassFragment[(int)'n']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'h']
            + g_staticParams.massUtility.pdAAMassFragment[(int)'h'];

   g_staticParams.massUtility.dCOminusH2 = g_staticParams.massUtility.dCO
            - g_staticParams.massUtility.pdAAMassFragment[(int)'h']
            - g_staticParams.massUtility.pdAAMassFragment[(int)'h'];

   GetParamValue("[COMET_ENZYME_INFO]", g_staticParams.enzymeInformation);
   if (!strncmp(g_staticParams.enzymeInformation.szSearchEnzymeBreakAA, "-", 1)
         && !strncmp(g_staticParams.enzymeInformation.szSearchEnzymeNoBreakAA, "-", 1))
   {
      g_staticParams.options.bNoEnzymeSelected = 1;
   }
   else
   {
      g_staticParams.options.bNoEnzymeSelected = 0;
   }

   GetParamValue("allowed_missed_cleavage", g_staticParams.enzymeInformation.iAllowedMissedCleavage);
   if (g_staticParams.enzymeInformation.iAllowedMissedCleavage < 0)
      g_staticParams.enzymeInformation.iAllowedMissedCleavage = 0;

   // Load ion series to consider, useA, useB, useY are for neutral losses.
   g_staticParams.ionInformation.iNumIonSeriesUsed = 0;
   for (int i=0; i<6; i++)
   {
      if (g_staticParams.ionInformation.iIonVal[i] > 0)
         g_staticParams.ionInformation.piSelectedIonSeries[g_staticParams.ionInformation.iNumIonSeriesUsed++] = i;
   }

   // Variable mod search for AAs listed in szVarModChar.
   g_staticParams.szMod[0] = '\0';
   g_staticParams.variableModParameters.bVarModSearch = false;
   g_staticParams.variableModParameters.bBinaryModSearch = false;
   g_staticParams.variableModParameters.bRequireVarMod = false;

   if (g_staticParams.peffInfo.iPeffSearch)
      g_staticParams.variableModParameters.bVarModSearch = true;

   for (int i=0; i<VMODS; i++)
   {
      if (!isEqual(g_staticParams.variableModParameters.varModList[i].dVarModMass, 0.0)
            && (g_staticParams.variableModParameters.varModList[i].szVarModChar[0]!='\0'))
      {
         sprintf(g_staticParams.szMod + strlen(g_staticParams.szMod), "(%s%c %+0.6f) ",
               g_staticParams.variableModParameters.varModList[i].szVarModChar,
               g_staticParams.variableModParameters.cModCode[i],
               g_staticParams.variableModParameters.varModList[i].dVarModMass);

         g_staticParams.variableModParameters.bVarModSearch = true;

         if (g_staticParams.variableModParameters.varModList[i].iBinaryMod)
            g_staticParams.variableModParameters.bBinaryModSearch = true;

         if (g_staticParams.variableModParameters.varModList[i].bRequireThisMod)
            g_staticParams.variableModParameters.bRequireVarMod = true;
      }
   }

   // Do Sp scoring after search based on how many lines to print out.
   if (g_staticParams.options.iNumStored < 1)
      g_staticParams.options.iNumStored = 1;

   if (g_staticParams.options.iNumPeptideOutputLines > g_staticParams.options.iNumStored)
      g_staticParams.options.iNumPeptideOutputLines = g_staticParams.options.iNumStored;
   else if (g_staticParams.options.iNumPeptideOutputLines < 1)
      g_staticParams.options.iNumPeptideOutputLines = 1;

   if (g_staticParams.peaksInformation.iNumMatchPeaks > 5)
      g_staticParams.peaksInformation.iNumMatchPeaks = 5;

   if (!isEqual(g_staticParams.staticModifications.dAddCterminusPeptide, 0.0))
   {
      sprintf(g_staticParams.szMod + strlen(g_staticParams.szMod), "+ct=%0.6f ",
            g_staticParams.staticModifications.dAddCterminusPeptide);
   }

   if (!isEqual(g_staticParams.staticModifications.dAddNterminusPeptide, 0.0))
   {
      sprintf(g_staticParams.szMod + strlen(g_staticParams.szMod), "+nt=%0.6f ",
            g_staticParams.staticModifications.dAddNterminusPeptide);
   }

   if (!isEqual(g_staticParams.staticModifications.dAddCterminusProtein, 0.0))
   {
      sprintf(g_staticParams.szMod + strlen(g_staticParams.szMod), "+ctprot=%0.6f ",
            g_staticParams.staticModifications.dAddCterminusProtein);
   }

   if (!isEqual(g_staticParams.staticModifications.dAddNterminusProtein, 0.0))
   {
      sprintf(g_staticParams.szMod + strlen(g_staticParams.szMod), "+ntprot=%0.6f ",
            g_staticParams.staticModifications.dAddNterminusProtein);
   }

   for (int i=65; i<=90; i++)  // 65-90 represents upper case letters in ASCII
   {
      if (!isEqual(g_staticParams.staticModifications.pdStaticMods[i], 0.0))
      {
         sprintf(g_staticParams.szMod + strlen(g_staticParams.szMod), "%c=%0.6f ", i,
               g_staticParams.massUtility.pdAAMassParent[i] += g_staticParams.staticModifications.pdStaticMods[i]);
         g_staticParams.massUtility.pdAAMassFragment[i] += g_staticParams.staticModifications.pdStaticMods[i];
      }
      else if (i=='B' || i=='J' || i=='X' || i=='Z')
      {
         g_staticParams.massUtility.pdAAMassParent[i] = 999999.;
         g_staticParams.massUtility.pdAAMassFragment[i] = 999999.;
      }
   }

   // Print out enzyme name to g_staticParams.szMod.
   if (!g_staticParams.options.bNoEnzymeSelected)
   {
      char szTmp[4];

      szTmp[0]='\0';
      if (g_staticParams.options.iEnzymeTermini != 2)
         sprintf(szTmp, ":%d", g_staticParams.options.iEnzymeTermini);

      sprintf(g_staticParams.szMod + strlen(g_staticParams.szMod), "Enzyme:%s (%d%s)",
            g_staticParams.enzymeInformation.szSearchEnzymeName,
            g_staticParams.enzymeInformation.iAllowedMissedCleavage,
            szTmp);
   }
   else
   {
      sprintf(g_staticParams.szMod + strlen(g_staticParams.szMod), "Enzyme:%s",
            g_staticParams.enzymeInformation.szSearchEnzymeName);
   }

   if (g_staticParams.tolerances.dFragmentBinStartOffset < 0.0
         || g_staticParams.tolerances.dFragmentBinStartOffset >1.0)
   {
      char szErrorMsg[256];
      sprintf(szErrorMsg,  " Error - bin offset %f must between 0.0 and 1.0\n",
            g_staticParams.tolerances.dFragmentBinStartOffset);
      string strErrorMsg(szErrorMsg);
      g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      return false;
   }

   if (!g_staticParams.options.bOutputOutFiles)
      g_staticParams.options.bSkipAlreadyDone = 0;

   g_staticParams.precalcMasses.dNtermProton = g_staticParams.staticModifications.dAddNterminusPeptide
      + PROTON_MASS;

   g_staticParams.precalcMasses.dCtermOH2Proton = g_staticParams.staticModifications.dAddCterminusPeptide
      + g_staticParams.massUtility.dOH2fragment
      + PROTON_MASS;

   g_staticParams.precalcMasses.dOH2ProtonCtermNterm = g_staticParams.massUtility.dOH2parent
      + PROTON_MASS
      + g_staticParams.staticModifications.dAddCterminusPeptide
      + g_staticParams.staticModifications.dAddNterminusPeptide;

   g_staticParams.enzymeInformation.iOneMinusOffset = 1 - g_staticParams.enzymeInformation.iSearchEnzymeOffSet;
   g_staticParams.enzymeInformation.iTwoMinusOffset = 2 - g_staticParams.enzymeInformation.iSearchEnzymeOffSet;

   string cutoff;
   FragmentCalculator& calculator = FragmentCalculator::instance();
   if (GetParamValue("relative_minimum_intensity", cutoff)) {
      calculator.setOption("relative_minimum_intensity", cutoff);
   }
   else {
       calculator.setOption("relative_minimum_intensity", "0.05");
   }

   string obsChargeOffset;
   if (GetParamValue("obs_charge_offset", obsChargeOffset)) {
      calculator.setOption("obs_charge_offset", obsChargeOffset);
   }

   string model_min_cutoff;
   if (GetParamValue("frag_isotope_min_prob", model_min_cutoff)) {
       calculator.setOption("frag_isotope_min_prob", model_min_cutoff);
   }

   string doNormalization;
   if (GetParamValue("normalization", doNormalization)) {
      calculator.setOption("normalization", doNormalization);
   }

   string doSqrt;
   if (GetParamValue("model_sqrt", doSqrt)) {
      calculator.setOption("model_sqrt", doSqrt);
   }

   string doWindowedScaling;
   if (GetParamValue("model_windowed_scaling", doWindowedScaling)) {
      calculator.setOption("model_windowed_scaling", doWindowedScaling);
   }

   string isUseIsotopeSqrtNormalization;
   if (GetParamValue("isotope_sqrt", isUseIsotopeSqrtNormalization)) {
       calculator.setOption("isotope_sqrt", isUseIsotopeSqrtNormalization);
   }

   string useMzRange;
   if (GetParamValue("use_mz_range", useMzRange)) {
      calculator.setOption("use_mz_range", useMzRange);
   }

   string modelType;
   if (GetParamValue("fragment_model", modelType)) {
      calculator.setOption("model", modelType);
   }
   else {
       calculator.setOption("model", "");
   }

   string isoMinProb;
   if (GetParamValue("isotope_min_prob", isoMinProb)) {
       calculator.setOption("frag_isotope_min_prob", isoMinProb);
   }

   string isotopeMassDiff;
   if (GetParamValue("isotope_mass_diff", isotopeMassDiff)) {
       calculator.setOption("isotope_mass_diff", isotopeMassDiff);
   }

   // Dist can only be set after model is selected.
   for (auto c : calculator.validOptionChars()) {
       string key = "dist-" + c;
       string value;
       if (GetParamValue(key, value)) {
           calculator.setOption(key, value);
       }
   }
   calculator.printOptions();
   return true;
}

void CometSearchManager::AddInputFiles(vector<InputFileInfo*> &pvInputFiles)
{
   int numInputFiles = (int)pvInputFiles.size();

   for (int i = 0; i < numInputFiles; i++)
      g_pvInputFiles.push_back(pvInputFiles.at(i));
}

void CometSearchManager::SetOutputFileBaseName(const char *pszBaseName)
{
   strcpy(g_staticParams.inputFile.szBaseName, pszBaseName);
}

std::map<std::string, CometParam*>& CometSearchManager::GetParamsMap()
{
   return _mapStaticParams;
}

void CometSearchManager::SetParam(const string &name, const string &strValue, const string& value)
{
   CometParam *pParam = new TypedCometParam<string>(CometParamType_String, strValue, value);
   pair<map<string, CometParam*>::iterator,bool> ret = _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   if (false == ret.second)
   {
      _mapStaticParams.erase(name);
      _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   }
}

bool CometSearchManager::GetParamValue(const string &name, string& value)
{
   std::map<string, CometParam*>::iterator it;
   it = _mapStaticParams.find(name);
   if (it == _mapStaticParams.end())
      return false;

   TypedCometParam<string> *pParam = static_cast<TypedCometParam<string>*>(it->second);
   value = pParam->GetValue();
   return true;
}

void CometSearchManager::SetParam(const std::string &name, const string &strValue, const int &value)
{
   CometParam *pParam = new TypedCometParam<int>(CometParamType_Int, strValue, value);
   pair<map<string, CometParam*>::iterator,bool> ret = _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   if (false == ret.second)
   {
      _mapStaticParams.erase(name);
      _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   }
}

bool CometSearchManager::GetParamValue(const string &name, int& value)
{
   std::map<string, CometParam*>::iterator it;
   it = _mapStaticParams.find(name);
   if (it == _mapStaticParams.end())
      return false;

   TypedCometParam<int> *pParam = static_cast<TypedCometParam<int>*>(it->second);
   value = pParam->GetValue();
   return true;
}

void CometSearchManager::SetParam(const std::string &name, const string &strValue, const long &value)
{
   CometParam *pParam = new TypedCometParam<long>(CometParamType_Long, strValue, value);
   pair<map<string, CometParam*>::iterator,bool> ret = _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   if (false == ret.second)
   {
      _mapStaticParams.erase(name);
      _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   }
}

bool CometSearchManager::GetParamValue(const string &name, long& value)
{
   std::map<string, CometParam*>::iterator it;
   it = _mapStaticParams.find(name);
   if (it == _mapStaticParams.end())
      return false;

   TypedCometParam<long> *pParam = static_cast<TypedCometParam<long>*>(it->second);
   value = pParam->GetValue();
   return true;
}

void CometSearchManager::SetParam(const string &name, const string &strValue, const double &value)
{
   CometParam *pParam = new TypedCometParam<double>(CometParamType_Double, strValue, value);
   pair<map<string, CometParam*>::iterator,bool> ret = _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   if (false == ret.second)
   {
      _mapStaticParams.erase(name);
      _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   }
}

bool CometSearchManager::GetParamValue(const string &name, double& value)
{
   std::map<string, CometParam*>::iterator it;
   it = _mapStaticParams.find(name);
   if (it == _mapStaticParams.end())
      return false;

   TypedCometParam<double> *pParam = static_cast<TypedCometParam<double>*>(it->second);
   value = pParam->GetValue();
   return true;
}

void CometSearchManager::SetParam(const string &name, const string &strValue, const VarMods &value)
{
   CometParam *pParam = new TypedCometParam<VarMods>(CometParamType_VarMods, strValue, value);
   pair<map<string, CometParam*>::iterator,bool> ret = _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   if (false == ret.second)
   {
      _mapStaticParams.erase(name);
      _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   }
}

bool CometSearchManager::GetParamValue(const string &name, VarMods & value)
{
   std::map<string, CometParam*>::iterator it;
   it = _mapStaticParams.find(name);
   if (it == _mapStaticParams.end())
      return false;

   TypedCometParam<VarMods> *pParam = static_cast<TypedCometParam<VarMods>*>(it->second);
   value = pParam->GetValue();
   return true;
}

void CometSearchManager::SetParam(const string &name, const string &strValue, const DoubleRange &value)
{
   CometParam *pParam = new TypedCometParam<DoubleRange>(CometParamType_DoubleRange, strValue, value);
   pair<map<string, CometParam*>::iterator,bool> ret = _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   if (false == ret.second)
   {
      _mapStaticParams.erase(name);
      _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   }
}

bool CometSearchManager::GetParamValue(const string &name, DoubleRange &value)
{
   std::map<string, CometParam*>::iterator it;
   it = _mapStaticParams.find(name);
   if (it == _mapStaticParams.end())
      return false;

   TypedCometParam<DoubleRange> *pParam = static_cast<TypedCometParam<DoubleRange>*>(it->second);
   value = pParam->GetValue();
   return true;
}

void CometSearchManager::SetParam(const string &name, const string &strValue, const IntRange &value)
{
   CometParam *pParam = new TypedCometParam<IntRange>(CometParamType_IntRange, strValue, value);
   pair<map<string, CometParam*>::iterator,bool> ret = _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   if (false == ret.second)
   {
      _mapStaticParams.erase(name);
      _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   }
}

bool CometSearchManager::GetParamValue(const string &name, IntRange &value)
{
   std::map<string, CometParam*>::iterator it;
   it = _mapStaticParams.find(name);
   if (it == _mapStaticParams.end())
      return false;

   TypedCometParam<IntRange> *pParam = static_cast<TypedCometParam<IntRange>*>(it->second);
   value = pParam->GetValue();
   return true;
}

void CometSearchManager::SetParam(const string &name, const string &strValue, const EnzymeInfo &value)
{
   CometParam *pParam = new TypedCometParam<EnzymeInfo>(CometParamType_EnzymeInfo, strValue, value);
   pair<map<string, CometParam*>::iterator,bool> ret = _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   if (false == ret.second)
   {
      _mapStaticParams.erase(name);
      _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   }
}

bool CometSearchManager::GetParamValue(const string &name, EnzymeInfo &value)
{
   std::map<string, CometParam*>::iterator it;
   it = _mapStaticParams.find(name);
   if (it == _mapStaticParams.end())
      return false;

   TypedCometParam<EnzymeInfo> *pParam = static_cast<TypedCometParam<EnzymeInfo>*>(it->second);
   value = pParam->GetValue();
   return true;
}

void CometSearchManager::SetParam(const string &name, const string &strValue, const vector<double> &value)
{
   CometParam *pParam = new TypedCometParam< vector<double> >(CometParamType_DoubleVector, strValue, value);
   pair<map<string, CometParam*>::iterator,bool> ret = _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   if (false == ret.second)
   {
      _mapStaticParams.erase(name);
      _mapStaticParams.insert(std::pair<std::string, CometParam*>(name, pParam));
   }
}

bool CometSearchManager::GetParamValue(const string &name,  vector<double> &value)
{
   std::map<string, CometParam*>::iterator it;
   it = _mapStaticParams.find(name);
   if (it == _mapStaticParams.end())
      return false;

   TypedCometParam< vector<double> > *pParam = static_cast<TypedCometParam< vector<double> >*>(it->second);
   value = pParam->GetValue();

   return true;
}

bool CometSearchManager::IsSearchError()
{
    return g_cometStatus.IsError();
}

void CometSearchManager::GetStatusMessage(string &strStatusMsg)
{
   g_cometStatus.GetStatusMsg(strStatusMsg);
}

bool CometSearchManager::IsValidCometVersion(const string &version)
{
    // Major version number must match to current binary
    if (strstr(comet_version, version.c_str()) || strstr("2018.01", version.c_str()) || strstr("2017.01", version.c_str())) //FIX
       return true;
    else
       return false;
}

void CometSearchManager::CancelSearch()
{
    g_cometStatus.SetStatus(CometResult_Cancelled, string("Search was cancelled."));
}

bool CometSearchManager::IsCancelSearch()
{
    return g_cometStatus.IsCancel();
}

void CometSearchManager::ResetSearchStatus()
{
    g_cometStatus.ResetStatus();
}

bool CometSearchManager::CreateIndex()
{
    // Override the Create Index flag to force it to create
    g_staticParams.options.bCreateIndex = 1;

    // The DoSearch will create the index and exit
    return DoSearch();
}

bool CometSearchManager::DoSearch()
{
   char szOut[256];

   if (!InitializeStaticParams())
      return false;

   if (!ValidateOutputFormat())
      return false;

   if (!ValidateSequenceDatabaseFile())
      return false;

   if (!ValidateScanRange())
      return false;

   if (!ValidatePeptideLengthRange())
      return false;

   bool bSucceeded = true;

   if (!g_staticParams.options.bOutputSqtStream && !g_staticParams.bIndexDb)
   {
      sprintf(szOut, " Comet version \"%s\"", comet_version);
//    if (!g_staticParams.options.bSkipUpdateCheck)
//       CometCheckForUpdates::CheckForUpdates(szOut);
      sprintf(szOut+strlen(szOut), "\n\n");

      logout(szOut);
      fflush(stdout);
   }

   if (g_staticParams.options.bCreateIndex) //index
   {
      bSucceeded = WriteIndexedDatabase();
      return bSucceeded;
   }

   if (g_staticParams.options.bOutputOutFiles)
      PrintOutfileHeader();

   bool bBlankSearchFile = false;

   if (strlen(g_staticParams.szDIAWindowsFile) > 0)
   {
      FILE *fp;

      if ((fp=fopen(g_staticParams.szDIAWindowsFile, "r")) != NULL)
      {
         // read DIA windows
         double dStartMass=0.0,
                dEndMass=0.0;
         char szTmp[512];

         while (fgets(szTmp, 512, fp))
         {
            sscanf(szTmp, "%lf %lf", &dStartMass, &dEndMass);
            if (dEndMass <= dStartMass)
            {
               char szErrorMsg[256];
               sprintf(szErrorMsg,  " Error - DIA window file end mass <= start mass:  %f %f.\n",  dStartMass, dEndMass);
               string strErrorMsg(szErrorMsg);
               g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
               logerr(szErrorMsg);

               g_pvDIAWindows.clear();
               return false;
            }

            g_pvDIAWindows.push_back(dStartMass);
            g_pvDIAWindows.push_back(dEndMass);
         }
         fclose(fp);
      }
      else
      {
         char szErrorMsg[256];
         sprintf(szErrorMsg,  " Error - cannot read DIA window file \"%s\".\n",  g_staticParams.szDIAWindowsFile);
         string strErrorMsg(szErrorMsg);
         g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
         logerr(szErrorMsg);
      }
   }

   for (int i=0; i<(int)g_pvInputFiles.size(); i++)
   {
      bSucceeded = UpdateInputFile(g_pvInputFiles.at(i));
      if (!bSucceeded)
         break;

      time_t tStartTime;
      time(&tStartTime);
      strftime(g_staticParams.szDate, 26, "%m/%d/%Y, %I:%M:%S %p", localtime(&tStartTime));

      if (!g_staticParams.options.bOutputSqtStream && !g_staticParams.bIndexDb)
      {
         sprintf(szOut, " Search start:  %s\n", g_staticParams.szDate);
         sprintf(szOut+strlen(szOut), " - Input file: %s\n", g_staticParams.inputFile.szFileName);
         logout(szOut);
         fflush(stdout);
      }

      int iFirstScan = g_staticParams.inputFile.iFirstScan;             // First scan to search specified by user.
      int iLastScan = g_staticParams.inputFile.iLastScan;               // Last scan to search specified by user.
      int iPercentStart = 0;                                            // percentage within input file for start scan of batch
      int iPercentEnd = 0;                                              // percentage within input file for end scan of batch
      int iAnalysisType = g_staticParams.inputFile.iAnalysisType;       // 1=dta (retired),
                                                                        // 2=specific scan,
                                                                        // 3=specific scan + charge,
                                                                        // 4=scan range,
                                                                        // 5=entire file

      // For SQT & pepXML output file, check if they can be written to before doing anything else.
      FILE *fpout_sqt=NULL;
      FILE *fpoutd_sqt=NULL;
      FILE *fpout_pepxml=NULL;
      FILE *fpoutd_pepxml=NULL;
      FILE *fpout_percolator=NULL;
      FILE *fpout_txt=NULL;
      FILE *fpoutd_txt=NULL;

      char szOutputSQT[1024];
      char szOutputDecoySQT[1024];
      char szOutputPepXML[1024];
      char szOutputDecoyPepXML[1024];
      char szOutputPercolator[1024];
      char szOutputTxt[1024];
      char szOutputDecoyTxt[1024];

      if (g_staticParams.options.bOutputSqtFile)
      {
         if (iAnalysisType == AnalysisType_EntireFile)
         {
#ifdef CRUX
            sprintf(szOutputSQT, "%s%s.target.sqt",
                  g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix);
#else
            sprintf(szOutputSQT, "%s%s.sqt",
                  g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix);
#endif
         }
         else
         {
#ifdef CRUX
            sprintf(szOutputSQT, "%s%s.%d-%d.target.sqt",
                  g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix, iFirstScan, iLastScan);
#else
            sprintf(szOutputSQT, "%s%s.%d-%d.sqt",
                  g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix, iFirstScan, iLastScan);
#endif
         }

         if ((fpout_sqt = fopen(szOutputSQT, "w")) == NULL)
         {
            char szErrorMsg[256];
            sprintf(szErrorMsg,  " Error - cannot write to file \"%s\".\n",  szOutputSQT);
            string strErrorMsg(szErrorMsg);
            g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
            logerr(szErrorMsg);
            bSucceeded = false;
         }

         CometWriteSqt::PrintSqtHeader(fpout_sqt, *this);

         if (bSucceeded && (g_staticParams.options.iDecoySearch == 2))
         {
            if (iAnalysisType == AnalysisType_EntireFile)
            {
               sprintf(szOutputDecoySQT, "%s%s.decoy.sqt",
                     g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix);
            }
            else
            {
               sprintf(szOutputDecoySQT, "%s%s.%d-%d.decoy.sqt",
                     g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix, iFirstScan, iLastScan);
            }

            if ((fpoutd_sqt = fopen(szOutputDecoySQT, "w")) == NULL)
            {
               char szErrorMsg[256];
               sprintf(szErrorMsg,  " Error - cannot write to decoy file \"%s\".\n",  szOutputDecoySQT);
               string strErrorMsg(szErrorMsg);
               g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
               logerr(szErrorMsg);
               bSucceeded = false;
            }

            CometWriteSqt::PrintSqtHeader(fpoutd_sqt, *this);
         }
      }

      if (bSucceeded && g_staticParams.options.bOutputTxtFile)
      {
         if (iAnalysisType == AnalysisType_EntireFile)
         {
#ifdef CRUX
            sprintf(szOutputTxt, "%s%s.target.txt",
                  g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix);
#else
            sprintf(szOutputTxt, "%s%s.txt",
                  g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix);
#endif
         }
         else
         {
#ifdef CRUX
            sprintf(szOutputTxt, "%s%s.%d-%d.target.txt",
                  g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix, iFirstScan, iLastScan);
#else
            sprintf(szOutputTxt, "%s%s.%d-%d.txt",
                  g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix, iFirstScan, iLastScan);
#endif
         }

         if ((fpout_txt = fopen(szOutputTxt, "w")) == NULL)
         {
            char szErrorMsg[256];
            sprintf(szErrorMsg,  " Error - cannot write to file \"%s\".\n",  szOutputTxt);
            string strErrorMsg(szErrorMsg);
            g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
            logerr(szErrorMsg);
            bSucceeded = false;
         }

         CometWriteTxt::PrintTxtHeader(fpout_txt);
         fflush(fpout_txt);

         if (bSucceeded && (g_staticParams.options.iDecoySearch == 2))
         {
            if (iAnalysisType == AnalysisType_EntireFile)
            {
               sprintf(szOutputDecoyTxt, "%s%s.decoy.txt",
                     g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix);
            }
            else
            {
               sprintf(szOutputDecoyTxt, "%s%s.%d-%d.decoy.txt",
                     g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix, iFirstScan, iLastScan);
            }

            if ((fpoutd_txt= fopen(szOutputDecoyTxt, "w")) == NULL)
            {
               char szErrorMsg[256];
               sprintf(szErrorMsg,  " Error - cannot write to decoy file \"%s\".\n",  szOutputDecoyTxt);
               string strErrorMsg(szErrorMsg);
               g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
               logerr(szErrorMsg);
               bSucceeded = false;
            }

            CometWriteTxt::PrintTxtHeader(fpoutd_txt);
         }
      }

      if (bSucceeded && g_staticParams.options.bOutputPepXMLFile)
      {
         if (iAnalysisType == AnalysisType_EntireFile)
         {
#ifdef CRUX
            sprintf(szOutputPepXML, "%s%s.target.pep.xml",
                  g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix);
#else
            sprintf(szOutputPepXML, "%s%s.pep.xml",
                  g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix);
#endif
         }
         else
         {
#ifdef CRUX
            sprintf(szOutputPepXML, "%s%s.%d-%d.target.pep.xml",
                  g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix, iFirstScan, iLastScan);
#else
            sprintf(szOutputPepXML, "%s%s.%d-%d.pep.xml",
                  g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix, iFirstScan, iLastScan);
#endif
         }

         if ((fpout_pepxml = fopen(szOutputPepXML, "w")) == NULL)
         {
            char szErrorMsg[256];
            sprintf(szErrorMsg,  " Error - cannot write to file \"%s\".\n",  szOutputPepXML);
            string strErrorMsg(szErrorMsg);
            g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
            logerr(szErrorMsg);
            bSucceeded = false;
         }

         if (bSucceeded)
            bSucceeded = CometWritePepXML::WritePepXMLHeader(fpout_pepxml, *this);

         if (bSucceeded && (g_staticParams.options.iDecoySearch == 2))
         {
            if (iAnalysisType == AnalysisType_EntireFile)
            {
               sprintf(szOutputDecoyPepXML, "%s%s.decoy.pep.xml",
                     g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix);
            }
            else
            {
               sprintf(szOutputDecoyPepXML, "%s%s.%d-%d.decoy.pep.xml",
                     g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix, iFirstScan, iLastScan);
            }

            if ((fpoutd_pepxml = fopen(szOutputDecoyPepXML, "w")) == NULL)
            {
               char szErrorMsg[256];
               sprintf(szErrorMsg,  " Error - cannot write to decoy file \"%s\".\n",  szOutputDecoyPepXML);
               string strErrorMsg(szErrorMsg);
               g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
               logerr(szErrorMsg);
               bSucceeded = false;
            }

            if (bSucceeded)
               bSucceeded = CometWritePepXML::WritePepXMLHeader(fpoutd_pepxml, *this);
         }
      }

      if (bSucceeded && g_staticParams.options.bOutputPercolatorFile)
      {
         if (iAnalysisType == AnalysisType_EntireFile)
         {
            sprintf(szOutputPercolator, "%s%s.pin",
                  g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix);
         }
         else
         {
            sprintf(szOutputPercolator, "%s%s.%d-%d.pin",
                  g_staticParams.inputFile.szBaseName, g_staticParams.szOutputSuffix, iFirstScan, iLastScan);
         }

         if ((fpout_percolator = fopen(szOutputPercolator, "w")) == NULL)
         {
            char szErrorMsg[256];
            sprintf(szErrorMsg,  " Error - cannot write to file \"%s\".\n",  szOutputPercolator);
            string strErrorMsg(szErrorMsg);
            g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
            logerr(szErrorMsg);
            bSucceeded = false;
         }

         if (bSucceeded)
         {
            // We need knowledge of max charge state in all searches
            // here in order to write the featureDescription header

            CometWritePercolator::WritePercolatorHeader(fpout_percolator);
         }
      }

      int iTotalSpectraSearched = 0;
      if (bSucceeded)
      {
         //MH: Allocate memory shared by threads during spectral processing.
         bSucceeded = CometPreprocess::AllocateMemory(g_staticParams.options.iNumThreads);
         if (!bSucceeded)
            break;

         // Allocate memory shared by threads during search
         bSucceeded = CometSearch::AllocateMemory(g_staticParams.options.iNumThreads);
         if (!bSucceeded)
            break;

         // For file access using MSToolkit.
         MSReader mstReader;

         // We want to read only MS2/MS3 scans.
         SetMSLevelFilter(mstReader);

         // We need to reset some of the static variables in-between input files
         CometPreprocess::Reset();

         int iBatchNum = 0;
         while (!CometPreprocess::DoneProcessingAllSpectra()) // Loop through iMaxSpectraPerSearch
         {
            iBatchNum++;
#ifdef PERF_DEBUG
            time_t tTotalSearchStartTime;
            time_t tTotalSearchEndTime;
            time_t tLoadAndPreprocessSpectraStartTime;
            time_t tLoadAndPreprocessSpectraEndTime;
            time_t tRunSearchStartTime;
            time_t tRunSearchEndTime;
            time_t tPostAnalysisStartTime;
            time_t tPostAnalysisEndTime;

            char szTimeBuffer[32];
            szTimeBuffer[0] = '\0';
#endif

            // Load and preprocess all the spectra.
            if (!g_staticParams.options.bOutputSqtStream && !g_staticParams.bIndexDb)
            {
               logout("   - Load spectra:");

#ifdef PERF_DEBUG
               char szOut[128];
               time(&tLoadAndPreprocessSpectraStartTime);
               strftime(szTimeBuffer, 26, "%m/%d/%Y, %I:%M:%S %p", localtime(&tLoadAndPreprocessSpectraStartTime));
               sprintf(szOut, " - Start LoadAndPreprocessSpectra:  %s\n", szTimeBuffer);
               logout(szOut);
#endif

               fflush(stdout);
            }

            g_cometStatus.SetStatusMsg(string("Loading and processing input spectra"));

            // IMPORTANT: From this point onwards, because we've loaded some
            // spectra, we MUST "goto cleanup_results" before exiting the loop,
            // or we will create a memory leak!

            bSucceeded = CometPreprocess::LoadAndPreprocessSpectra(mstReader,
                iFirstScan, iLastScan, iAnalysisType,
                g_staticParams.options.iNumThreads,  // min # threads
                g_staticParams.options.iNumThreads); // max # threads

            if (!bSucceeded)
               goto cleanup_results;

            iPercentStart = iPercentEnd;
            iPercentEnd = mstReader.getPercent();

#ifdef PERF_DEBUG
            if (!g_staticParams.options.bOutputSqtStream)
            {
               char szOut[128];
               time(&tLoadAndPreprocessSpectraEndTime);
               strftime(szTimeBuffer, 26, "%m/%d/%Y, %I:%M:%S %p", localtime(&tLoadAndPreprocessSpectraEndTime));
               sprintf(szOut, " - End LoadAndPreprocessSpectra:  %s\n", szTimeBuffer);
               logout(szOut);
               int iElapsedTime=(int)difftime(tLoadAndPreprocessSpectraEndTime, tLoadAndPreprocessSpectraStartTime);
               sprintf(szOut, " - Time spent in LoadAndPreprocessSpectra:  %d seconds\n", iElapsedTime);
               logout(szOut);
               fflush(stdout);
            }
#endif

            if (g_pvQuery.empty())
               continue;    //FIX make sure continue instead of break makes sense
                            // possible no spectrum in batch passes filters; do not want to break in that case;
 //              break; // no search to run
            else
               iTotalSpectraSearched += (int)g_pvQuery.size();

            bSucceeded = AllocateResultsMem();

            if (!bSucceeded)
               goto cleanup_results;

            char szStatusMsg[256];
            sprintf(szStatusMsg, " %d\n", (int)g_pvQuery.size());
            if (!g_staticParams.options.bOutputSqtStream && !g_staticParams.bIndexDb)
            {
               char szOut[128];
               sprintf(szOut, "%s", szStatusMsg);
               logout(szOut);
            }
            g_cometStatus.SetStatusMsg(string(szStatusMsg));

            if (g_staticParams.options.bMango)
            {
               int iCurrentScanNumber = 0;       // used to track multiple Mango precursors from same scan number
               int iMangoIndex=0;

               // sort back to original spectrum order in MS2 scan in order to associate pairs
               // based on sequential order of precursors for each scan
               std::sort(g_pvQuery.begin(), g_pvQuery.end(), compareByMangoIndex);

               for (std::vector<Query*>::iterator it = g_pvQuery.begin(); it != g_pvQuery.end(); ++it)
               {
                  if ((*it)->_spectrumInfoInternal.iScanNumber != iCurrentScanNumber)
                  {
                     iCurrentScanNumber = (*it)->_spectrumInfoInternal.iScanNumber;
                     iMangoIndex = 0;
                  }
                  else
                     iMangoIndex++;

                  sprintf((*it)->_spectrumInfoInternal.szMango, "%03d_%c", (int)iMangoIndex/2, (iMangoIndex % 2)?'B':'A');
               }
            }

            if (g_pvDIAWindows.size() > 0)
            {
               // delete scan entries that do not map to a DIA window;
               // should we throw error and abort search in this case instead?
               std::vector<Query*>::iterator it = g_pvQuery.begin();

               while(it != g_pvQuery.end())
               {
                  if ( (*it)->_pepMassInfo.dPeptideMassToleranceMinus == 0.0
                        || (*it)->_pepMassInfo.dPeptideMassTolerancePlus == 0.0)
                  {
                     g_pvQuery.erase(it);
                  }
                  else
                     ++it;
               }
            }

            // Sort g_pvQuery vector by dExpPepMass.
            std::sort(g_pvQuery.begin(), g_pvQuery.end(), compareByPeptideMass);

            g_massRange.dMinMass = g_pvQuery.at(0)->_pepMassInfo.dPeptideMassToleranceMinus;
            g_massRange.dMaxMass = g_pvQuery.at(g_pvQuery.size()-1)->_pepMassInfo.dPeptideMassTolerancePlus;
            if(g_massRange.dMaxMass - g_massRange.dMinMass > g_massRange.dMinMass)
               g_massRange.bNarrowMassRange = true;
            else
               g_massRange.bNarrowMassRange = false;

#ifdef PERF_DEBUG
            if (!g_staticParams.options.bOutputSqtStream)
            {
               char szOut[128];
               time(&tRunSearchStartTime);
               strftime(szTimeBuffer, 26, "%m/%d/%Y, %I:%M:%S %p", localtime(&tRunSearchStartTime));
               sprintf(szOut, " - Start RunSearch:  %s\n", szTimeBuffer);
               logout(szOut);
               fflush(stdout);
            }
#endif

            bSucceeded = !g_cometStatus.IsError() && !g_cometStatus.IsCancel();
            if (!bSucceeded)
               goto cleanup_results;

            g_cometStatus.SetStatusMsg(string("Running search..."));

            // Now that spectra are loaded to memory and sorted, do search.
            bSucceeded = CometSearch::RunSearch(g_staticParams.options.iNumThreads, g_staticParams.options.iNumThreads, iPercentStart, iPercentEnd);
            if (!bSucceeded)
               goto cleanup_results;

#ifdef PERF_DEBUG
            if (!g_staticParams.options.bOutputSqtStream)
            {
               char szOut[128];
               time(&tRunSearchEndTime);
               strftime(szTimeBuffer, 26, "%m/%d/%Y, %I:%M:%S %p", localtime(&tRunSearchEndTime));
               sprintf(szOut, " - End RunSearch:  %s\n", szTimeBuffer);
               logout(szOut);

               int iElapsedTime=(int)difftime(tRunSearchEndTime, tRunSearchStartTime);
               sprintf(szOut, " - Time spent in RunSearch:  %d seconds\n", iElapsedTime);
               logout(szOut);

               time(&tPostAnalysisStartTime);
               strftime(szTimeBuffer, 26, "%m/%d/%Y, %I:%M:%S %p", localtime(&tPostAnalysisStartTime));
               sprintf(szOut, " - Start PostAnalysis:  %s\n", szTimeBuffer);
               logout(szOut);

               fflush(stdout);
            }
#endif

            bSucceeded = !g_cometStatus.IsError() && !g_cometStatus.IsCancel();
            if (!bSucceeded)
               goto cleanup_results;

            if (!g_staticParams.options.bOutputSqtStream && !g_staticParams.bIndexDb)
            {
               logout("     - Post analysis:");
               fflush(stdout);
            }

            g_cometStatus.SetStatusMsg(string("Performing post-search analysis ..."));

            // Sort each entry by xcorr, calculate E-values, etc.
            bSucceeded = CometPostAnalysis::PostAnalysis(g_staticParams.options.iNumThreads, g_staticParams.options.iNumThreads);
            if (!bSucceeded)
               goto cleanup_results;

#ifdef PERF_DEBUG
            if (!g_staticParams.options.bOutputSqtStream)
            {
               char szOut[128];
               time(&tPostAnalysisEndTime);
               strftime(szTimeBuffer, 26, "%m/%d/%Y, %I:%M:%S %p", localtime(&tPostAnalysisEndTime));
               sprintf(szOut, " - End PostAnalysis:  %s\n", szTimeBuffer);
               logout(szOut);
               int iElapsedTime=(int)difftime(tPostAnalysisEndTime, tPostAnalysisStartTime);
               sprintf(szOut, " - Time spent in PostAnalysis:  %d seconds\n", iElapsedTime);
               logout(szOut);
               fflush(stdout);
            }
#endif

            // Sort g_pvQuery vector by scan.
            std::sort(g_pvQuery.begin(), g_pvQuery.end(), compareByScanNumber);

            CalcRunTime(tStartTime);

            FILE *fpdb;  // need FASTA file again to grab headers for output (currently just store file positions)
            if ((fpdb=fopen(g_staticParams.databaseInfo.szDatabase, "rb")) == NULL)
            {
               char szErrorMsg[1024];
               sprintf(szErrorMsg, " Error - cannot read database file \"%s\".\n", g_staticParams.databaseInfo.szDatabase);
               string strErrorMsg(szErrorMsg);
               g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
               logerr(szErrorMsg);
               return false;
            }

            if (!g_staticParams.options.bOutputSqtStream && !g_staticParams.bIndexDb)
            {
               logout("  done\n");
               fflush(stdout);
            }

            if (g_staticParams.options.bOutputOutFiles)
            {
               bSucceeded = CometWriteOut::WriteOut(fpdb);
               if (!bSucceeded)
                  goto cleanup_results;
            }

            if (g_staticParams.options.bOutputPepXMLFile)
               CometWritePepXML::WritePepXML(fpout_pepxml, fpoutd_pepxml, fpdb);

            if (g_staticParams.options.bOutputPercolatorFile)
            {
               bSucceeded = CometWritePercolator::WritePercolator(fpout_percolator, fpdb);
               if (!bSucceeded)
                  goto cleanup_results;
            }

            if (g_staticParams.options.bOutputTxtFile)
               CometWriteTxt::WriteTxt(fpout_txt, fpoutd_txt, fpdb);

            //// Write SQT last as I destroy the g_staticParams.szMod string during that process
            if (g_staticParams.options.bOutputSqtStream || g_staticParams.options.bOutputSqtFile)
               CometWriteSqt::WriteSqt(fpout_sqt, fpoutd_sqt, fpdb);

            fclose(fpdb);

   cleanup_results:
            // Deleting each Query object in the vector calls its destructor, which
            // frees the spectral memory (see definition for Query in CometData.h).
            for (std::vector<Query*>::iterator it = g_pvQuery.begin(); it != g_pvQuery.end(); ++it)
               delete *it;

            g_pvQuery.clear();

            if (!bSucceeded)
               break;
         }

         if (bSucceeded)
         {
            if (iTotalSpectraSearched == 0)
               logout(" Warning - no spectra searched.\n\n");

            if (!g_staticParams.options.bOutputSqtStream && !g_staticParams.bIndexDb)
            {
               char szOut[128];
               time_t tEndTime;

               time(&tEndTime);
               int iElapsedTime = (int)difftime(tEndTime, tStartTime);

               strftime(g_staticParams.szDate, 26, "%m/%d/%Y, %I:%M:%S %p", localtime(&tEndTime));
               sprintf(szOut, " Search end:    %s", g_staticParams.szDate);

               int hours, mins, secs;

               hours = (int)(iElapsedTime/3600);
               mins = (int)(iElapsedTime/60) - (hours*60);
               secs = (int)(iElapsedTime%60);

               if (hours)
                  sprintf(szOut+strlen(szOut), ", %dh:%dm:%ds", hours, mins, secs);
               else
                  sprintf(szOut+strlen(szOut), ", %dm:%ds", mins, secs);
               sprintf(szOut+strlen(szOut), "\n\n");

               logout(szOut);
            }

            if (NULL != fpout_pepxml)
               CometWritePepXML::WritePepXMLEndTags(fpout_pepxml);

            if (NULL != fpoutd_pepxml)
               CometWritePepXML::WritePepXMLEndTags(fpoutd_pepxml);
         }
      }

      // Clean up the input files vector
      g_staticParams.vectorMassOffsets.clear();

      //MH: Deallocate spectral processing memory.
      CometPreprocess::DeallocateMemory(g_staticParams.options.iNumThreads);

      // Deallocate search memory
      CometSearch::DeallocateMemory(g_staticParams.options.iNumThreads);

      if (NULL != fpout_pepxml)
      {
         fclose(fpout_pepxml);
         fpout_pepxml = NULL;
         if (iTotalSpectraSearched == 0)
            unlink(szOutputPepXML);
      }

      if (NULL != fpoutd_pepxml)
      {
         fclose(fpoutd_pepxml);
         fpoutd_pepxml = NULL;
         if (iTotalSpectraSearched == 0)
            unlink(szOutputDecoyPepXML);
      }

      if (NULL != fpout_percolator)
      {
         fclose(fpout_percolator);
         fpout_percolator = NULL;
         if (iTotalSpectraSearched == 0)
            unlink(szOutputPercolator);
      }

      if (NULL != fpout_sqt)
      {
         fclose(fpout_sqt);
         fpout_sqt = NULL;
         if (iTotalSpectraSearched == 0)
            unlink(szOutputSQT);
      }

      if (NULL != fpoutd_sqt)
      {
         fclose(fpoutd_sqt);
         fpoutd_sqt = NULL;
         if (iTotalSpectraSearched == 0)
            unlink(szOutputDecoySQT);
      }

      if (NULL != fpout_txt)
      {
         fclose(fpout_txt);
         fpout_txt = NULL;
         if (iTotalSpectraSearched == 0)
            unlink(szOutputTxt);
      }

      if (NULL != fpoutd_txt)
      {
         fclose(fpoutd_txt);
         fpoutd_txt = NULL;
         if (iTotalSpectraSearched == 0)
            unlink(szOutputDecoyTxt);
      }

      if (iTotalSpectraSearched == 0)
         bBlankSearchFile = true;

      if (!bSucceeded)
         break;
   }

   if (bBlankSearchFile)
      return true;
   else
      return bSucceeded;
}


bool CometSearchManager::InitializeSingleSpectrumSearch()
{
   // Skip doing if already completed successfully.
   if (singleSearchInitializationComplete)
   {
      return true;
   }

   if (!InitializeStaticParams())
      return false;

   if (!ValidateSequenceDatabaseFile())
       return false;

   // this uses a single thread
   singleSearchThreadCount = 1;
   g_staticParams.options.iNumThreads = singleSearchThreadCount;

   bool bSucceeded;
   //MH: Allocate memory shared by threads during spectral processing.
   bSucceeded = CometPreprocess::AllocateMemory(g_staticParams.options.iNumThreads);
   if (!bSucceeded)
       return bSucceeded;

    // Allocate memory shared by threads during search
   bSucceeded = CometSearch::AllocateMemory(g_staticParams.options.iNumThreads);
   if (!bSucceeded)
      return bSucceeded;

   singleSearchInitializationComplete = true;
   return true;
}


void CometSearchManager::FinalizeSingleSpectrumSearch()
{
   if (singleSearchInitializationComplete)
   {
      //MH: Deallocate spectral processing memory.
      CometPreprocess::DeallocateMemory(singleSearchThreadCount);

      // Deallocate search memory
      CometSearch::DeallocateMemory(singleSearchThreadCount);

      singleSearchInitializationComplete = false;
   }
}


int CometSearchManager::CheckIdxPrecursorMatch(const int iPrecursorCharge,
                                               const double dMZ)
{
   // Now that spectra are loaded to memory and sorted, do search.
   int iReturn;

   CometSearch::CheckIdxPrecursorMatch(iPrecursorCharge, dMZ, &iReturn);

   // if iReturn == 1, the mz/charge is found in the index database

   return iReturn;
}

bool CometSearchManager::DoSingleSpectrumSearch(int iPrecursorCharge,
                                                double dMZ,
                                                double* pdMass,
                                                double* pdInten,
                                                int iNumPeaks,
                                                string& strReturnPeptide,
                                                string& strReturnProtein,
                                                vector<Fragment> & matchedFragments,
                                                Scores & score)
{
   int iPercentStart = 0;
   int iPercentEnd = 0;

   if (iNumPeaks == 0)
      return false;

   if (dMZ * iPrecursorCharge - (iPrecursorCharge - 1)*PROTON_MASS > g_staticParams.options.dPeptideMassHigh)
   {
      // this assumes dPeptideMassHigh is set correctly in the calling program
      return false;
   }

   if (!InitializeSingleSpectrumSearch())
   {
      return false;
   }

   // We need to reset some of the static variables in-between input files
   CometPreprocess::Reset();

   // IMPORTANT: From this point onwards, because we've loaded some
   // spectra, we MUST "goto cleanup_results" before exiting the loop,
   // or we will create a memory leak!

   int iArraySize = (int)((g_staticParams.options.dPeptideMassHigh + g_staticParams.tolerances.dInputTolerance + 2.0) * g_staticParams.dInverseBinWidth);
   double *pdTmpSpectrum = new double[iArraySize];  // use this to determine most intense b/y-ions masses to report back
   bool bSucceeded = CometPreprocess::PreprocessSingleSpectrum(iPrecursorCharge, dMZ, pdMass, pdInten, iNumPeaks, pdTmpSpectrum);

   if (!bSucceeded)
      goto cleanup_results;

   if (g_pvQuery.empty())
   {
      delete[] pdTmpSpectrum;
      return false; // no search to run
   }

   bSucceeded = AllocateResultsMem();

   if (!bSucceeded)
      goto cleanup_results;

   g_massRange.dMinMass = g_pvQuery.at(0)->_pepMassInfo.dPeptideMassToleranceMinus;
   g_massRange.dMaxMass = g_pvQuery.at(g_pvQuery.size()-1)->_pepMassInfo.dPeptideMassTolerancePlus;

   if (g_massRange.dMaxMass - g_massRange.dMinMass > g_massRange.dMinMass)
      g_massRange.bNarrowMassRange = true;  // unused in this context but setting here anyways
   else
      g_massRange.bNarrowMassRange = false;

   // Now that spectra are loaded to memory and sorted, do search.
   bSucceeded = CometSearch::RunSearch(g_staticParams.options.iNumThreads, g_staticParams.options.iNumThreads, iPercentStart, iPercentEnd);

   if (bSucceeded && g_pvQuery.at(0)->iMatchPeptideCount > 0)
      CometPostAnalysis::AnalyzeSP(0);
   else
      goto cleanup_results;

   bSucceeded = !g_cometStatus.IsError() && !g_cometStatus.IsCancel();

   if (!bSucceeded)
      goto cleanup_results;

   int iSize;
   iSize  = g_pvQuery.at(0)->iMatchPeptideCount;
   if (iSize > g_staticParams.options.iNumStored)
      iSize = g_staticParams.options.iNumStored;

   // simply take top xcorr peptide as E-value calculation too expensive
   if (iSize > 1) {
      std::sort(g_pvQuery.at(0)->_pResults, g_pvQuery.at(0)->_pResults + iSize, CometPostAnalysis::SortFnXcorr);
   }

   Query* pQuery;
   pQuery = g_pvQuery.at(0);  // return info for top hit only

   if (iSize > 0 && pQuery->_pResults[0].fXcorr>0.0 && pQuery->_pResults[0].iLenPeptide>0)
   {
      Results *pOutput = pQuery->_pResults;

      // Set return values for peptide sequence, protein, xcorr and E-value
      strReturnPeptide = std::string(1, pOutput[0].szPrevNextAA[0]) + ".";
      for (int i=0; i< pOutput[0].iLenPeptide; i++)
      {
         strReturnPeptide += pOutput[0].szPeptide[i];

         if (pOutput[0].piVarModSites[i] != 0)
         {
            std::stringstream ss;
            ss << "[" << std::fixed << std::setprecision(4) << pOutput[0].pdVarModSites[i] << "]";
            strReturnPeptide += ss.str();
         }
      }
      strReturnPeptide += "." + std::string(1, pOutput[0].szPrevNextAA[1]);

      strReturnProtein = pOutput[0].strSingleSearchProtein;            //protein

      score.xCorr         = pOutput[0].fXcorr;                        // xcorr
      score.mass          = pOutput[0].dPepMass - PROTON_MASS;        // calc neutral pep mass
      score.matchedIons   = pOutput[0].iMatchedIons;                  // ions matched
      score.totalIons     = pOutput[0].iTotalIons;                    // ions tot

      int iMinLength = MAX_PEPTIDE_LEN;
      for (int x = 0; x < iSize; x++)
      {
         int iLen = (int)strlen(pOutput[x].szPeptide);
         if (iLen == 0)
            break;
         if (iLen < iMinLength)
            iMinLength = iLen;
      }

      double dDeltaCn = 1.0;       // this is deltaCn between top hit and peptide in list (or next dissimilar peptide)

      if (iSize > 1)
      {
         for (int j = 1; j < iSize; j++)
         {
            // very poor way of calculating peptide similarity but it's what we have for now
            int iDiffCt = 0;

            for (int k = 0; k < iMinLength; k++)
            {
               // I-L and Q-K are same for purposes here
               if (pOutput[0].szPeptide[k] != pOutput[j].szPeptide[k])
               {
                  if (!((pOutput[0].szPeptide[k] == 'K' || pOutput[0].szPeptide[k] == 'Q')
                     && (pOutput[j].szPeptide[k] == 'K' || pOutput[j].szPeptide[k] == 'Q'))
                     && !((pOutput[0].szPeptide[k] == 'I' || pOutput[0].szPeptide[k] == 'L')
                        && (pOutput[j].szPeptide[k] == 'I' || pOutput[j].szPeptide[k] == 'L')))
                  {
                     iDiffCt++;
                  }
               }
            }

            // calculate deltaCn only if sequences are less than 0.75 similar
            if (((double)(iMinLength - iDiffCt) / iMinLength) < 0.75)
            {
               if (pOutput[0].fXcorr > 0.0 && pOutput[j].fXcorr >= 0.0)
                  dDeltaCn = 1.0 - pOutput[j].fXcorr / pOutput[0].fXcorr;
               else if (pOutput[0].fXcorr > 0.0 && pOutput[j].fXcorr < 0.0)
                  dDeltaCn = 1.0;
               else
                  dDeltaCn = 0.0;
               break;
            }
         }
      }

      score.dCn = dDeltaCn;  // dCn is calculated from first dis-similar peptide

      // Conversion table from b/y ions to the other types (a,c,x,z)
      const double ionMassesRelative[NUM_ION_SERIES] =
      {
         // N term relative
         -(Carbon_Mono + Oxygen_Mono), // a (CO difference from b)
         0, // b
         (Nitrogen_Mono + (3 * Hydrogen_Mono)), // c (NH3 difference from b)

         // C Term relative
         (Carbon_Mono + Oxygen_Mono - (2 * Hydrogen_Mono)), // x (CO-2H difference from y)
          0, // y
         -(Nitrogen_Mono + (2 * Hydrogen_Mono)), // z (NH2 difference from y)

         // Not Used
         0, // not used
         0, // not used
         0  // not used
      };

      // now deal with calculating b- and y-ions and returning most intense matches
      double dBion = g_staticParams.precalcMasses.dNtermProton;
      double dYion = g_staticParams.precalcMasses.dCtermOH2Proton;

      if (pQuery->_pResults[0].szPrevNextAA[0] == '-')
      {
         dBion += g_staticParams.staticModifications.dAddNterminusProtein;
      }
      if (pQuery->_pResults[0].szPrevNextAA[1] == '-')
      {
         dYion += g_staticParams.staticModifications.dAddCterminusProtein;
      }

      // mods at peptide length +1 and +2 are for n- and c-terminus
      if (g_staticParams.variableModParameters.bVarModSearch
         && (pQuery->_pResults[0].piVarModSites[pQuery->_pResults[0].iLenPeptide] == 1))
      {
         dBion += g_staticParams.variableModParameters.varModList[pQuery->_pResults[0].piVarModSites[pQuery->_pResults[0].iLenPeptide] - 1].dVarModMass;
      }

      if (g_staticParams.variableModParameters.bVarModSearch
         && (pQuery->_pResults[0].piVarModSites[pQuery->_pResults[0].iLenPeptide + 1] == 1))
      {
         dYion += g_staticParams.variableModParameters.varModList[pQuery->_pResults[0].piVarModSites[pQuery->_pResults[0].iLenPeptide + 1] - 1].dVarModMass;
      }

      //vector<MatchedIonsStruct> vMatchedYions;
      //vector<MatchedIonsStruct> vMatchedBions;

      int iTmp;

      // Generate pdAAforward for pQuery->_pResults[0].szPeptide.
      for (int i = 0; i < pQuery->_pResults[0].iLenPeptide - 1; i++)
      {
         int iPos = pQuery->_pResults[0].iLenPeptide - i - 1;

         dBion += g_staticParams.massUtility.pdAAMassFragment[(int)pQuery->_pResults[0].szPeptide[i]];
         dYion += g_staticParams.massUtility.pdAAMassFragment[(int)pQuery->_pResults[0].szPeptide[iPos]];

         if (g_staticParams.variableModParameters.bVarModSearch)
         {
            if (pQuery->_pResults[0].piVarModSites[i] != 0)
               dBion += pQuery->_pResults[0].pdVarModSites[i];

            if (pQuery->_pResults[0].piVarModSites[iPos] != 0)
               dYion += pQuery->_pResults[0].pdVarModSites[iPos];
         }

         for (int ctCharge = 1; ctCharge <= pQuery->_spectrumInfoInternal.iMaxFragCharge; ctCharge++)
         {
            // calculate every ion series the user specified
            for (int ionSeries = 0; ionSeries <= ION_SERIES_Z; ionSeries++)
            {
               // skip ion series that are not enabled.
               if (!g_staticParams.ionInformation.iIonVal[ionSeries])
               {
                  continue;
               }

               bool isNTerm = (ionSeries <= ION_SERIES_C);

               // get the fragment mass if it is n- or c-terimnus
               double mass = (isNTerm) ? dBion : dYion;
               int fragNumber = i + 1;

               // Add any conversion factor from differnt ion series (e.g. b -> a, or y -> z)
               mass += ionMassesRelative[ionSeries];

               double mz = (mass + (ctCharge - 1)*PROTON_MASS) / ctCharge;

               iTmp = BIN(mz);
               if (iTmp<iArraySize && pdTmpSpectrum[iTmp] > 0.0)
               {
                  Fragment frag;
                  frag.intensity = pdTmpSpectrum[iTmp];
                  frag.mass = mass;
                  frag.type = ionSeries;
                  frag.number = fragNumber;
                  frag.charge = ctCharge;
                  matchedFragments.push_back(frag);
               }
            }
         }
      }
   }
   else
   {
      strReturnPeptide = "";  // peptide
      strReturnProtein = "";  // protein
      score.xCorr         = -1;       // xcorr
      score.mass          = 0;        // calc neutral pep mass
      score.matchedIons   = 0;        // ions matched
      score.totalIons     = 0;        // ions tot
      score.dCn           = 0;        // dCn
   }

cleanup_results:

   // Deleting each Query object in the vector calls its destructor, which
   // frees the spectral memory (see definition for Query in CometData.h).
   for (std::vector<Query*>::iterator it = g_pvQuery.begin(); it != g_pvQuery.end(); ++it)
      delete *it;

   g_pvQuery.clear();

   // Clean up the input files vector
   g_staticParams.vectorMassOffsets.clear();

   delete[] pdTmpSpectrum;

   return bSucceeded;
}


// this sort function needs to compare peptide, modification state, and protein file pointer
bool CometSearchManager::CompareByPeptide(const DBIndex &lhs,
                                          const DBIndex &rhs)
{
   if (!strcmp(lhs.szPeptide, rhs.szPeptide))
   {
      // peptides are same here so look at mass next
      if (fabs(lhs.dPepMass - rhs.dPepMass) > FLOAT_ZERO)
      {
         // masses are different
         if (lhs.dPepMass < rhs.dPepMass)
            return true;
         else
            return false;
      }

      // same sequences and masses here so next look at mod state
      for (unsigned int i=0; i<strlen(lhs.szPeptide)+2; i++)
      {
         if (lhs.pcVarModSites[i] != rhs.pcVarModSites[i])
         {
            // different mod state
            if (lhs.pcVarModSites[i] > rhs.pcVarModSites[i])
               return true;
            else
               return false;
         }
      }

      // at this point, same peptide, same mass, same mods so return first protein
      if (lhs.lIndexProteinFilePosition < rhs.lIndexProteinFilePosition)
         return true;
      else
         return false;
   }

   // peptides are different
   if (strcmp(lhs.szPeptide, rhs.szPeptide)<0)
      return true;
   else
      return false;
};


// sort by mass, then peptide, then modification state, then protein fp location
bool CometSearchManager::CompareByMass(const DBIndex &lhs,
                                       const DBIndex &rhs)
{
   if (fabs(lhs.dPepMass - rhs.dPepMass) > FLOAT_ZERO)
   {
      // masses are different
      if (lhs.dPepMass < rhs.dPepMass)
         return true;
      else
         return false;
   }

   // at this point, peptides are same mass so next need to compare sequences

   if (!strcmp(lhs.szPeptide, rhs.szPeptide))
   {
      // same sequences and masses here so next look at mod state
      for (unsigned int i=0; i<strlen(lhs.szPeptide)+2; i++)
      {
         if (lhs.pcVarModSites[i] != rhs.pcVarModSites[i])
         {
            if (lhs.pcVarModSites[i] > rhs.pcVarModSites[i])
               return true;
            else
               return false;
         }
      }

      // at this point, same peptide, same mass, same mods so return first protein
      if (lhs.lIndexProteinFilePosition < rhs.lIndexProteinFilePosition)
         return true;
      else
         return false;
   }

   // if here, peptide sequences are different (but w/same mass) so sort alphabetically
   if (strcmp(lhs.szPeptide, rhs.szPeptide) < 0)
      return true;
   else
      return false;

}


bool CometSearchManager::WriteIndexedDatabase(void)
{
   FILE *fptr;
   bool bSucceeded;
   char szOut[256];

   char szIndexFile[SIZE_FILE];
   sprintf(szIndexFile, "%s.idx", g_staticParams.databaseInfo.szDatabase);

   if ((fptr = fopen(szIndexFile, "wb")) == NULL)
   {
      printf(" Error - cannot open index file %s to write\n", szIndexFile);
      exit(1);
   }

   sprintf(szOut, " Creating peptide index file ");
   logout(szOut);
   fflush(stdout);

   bSucceeded = CometSearch::AllocateMemory(g_staticParams.options.iNumThreads);

   g_massRange.dMinMass = g_staticParams.options.dPeptideMassLow;
   g_massRange.dMaxMass = g_staticParams.options.dPeptideMassHigh;

   if (g_massRange.dMaxMass - g_massRange.dMinMass > g_massRange.dMinMass)
      g_massRange.bNarrowMassRange = true;
   else
      g_massRange.bNarrowMassRange = false;

   if (bSucceeded)
       bSucceeded = CometSearch::RunSearch(g_staticParams.options.iNumThreads, g_staticParams.options.iNumThreads, 0, 0);

   if (!bSucceeded)
   {
      char szErrorMsg[256];
      sprintf(szErrorMsg, " Error performing RunSearch() to create indexed database. \n");
      logerr(szErrorMsg);
      CometSearch::DeallocateMemory(g_staticParams.options.iNumThreads);
      return false;
   }

   // sanity check
   if (g_pvDBIndex.size() == 0)
   {
      char szErrorMsg[256];
      sprintf(szErrorMsg, " Error - no peptides in index; check the input database file.\n");
      logerr(szErrorMsg);
      CometSearch::DeallocateMemory(g_staticParams.options.iNumThreads);
      return false;
   }

   // remove duplicates
   sprintf(szOut, " - removing duplicates\n");
   logout(szOut);
   fflush(stdout);

   // keep unique entries only; sort by peptide/modification state and protein
   sort(g_pvDBIndex.begin(), g_pvDBIndex.end(), CompareByPeptide);

   // At this point, need to create g_pvProteinsList protein file position vector of vectors to map each peptide
   // to every protein. g_pvDBIndex.at().lProteinFilePosition is now reference to protein vector entry
   long lProtCount = 0;
   vector<vector<comet_fileoffset_t> > g_pvProteinsList;
   g_pvProteinsList.clear();
   for (size_t i = 0; i < g_pvDBIndex.size(); i++)
   {
      if (i == 0)
      {
         //store new protein offset containing first protein file position
         vector<comet_fileoffset_t> temp;
         temp.push_back(g_pvDBIndex.at(i).lIndexProteinFilePosition);
         g_pvProteinsList.push_back(temp);
         g_pvDBIndex.at(i).lIndexProteinFilePosition = lProtCount;
         lProtCount++;
      }
      else
      {
         if (!strcmp(g_pvDBIndex.at(i).szPeptide, g_pvDBIndex.at(i-1).szPeptide))
         {
            // peptide is same as previous, add in the current protein offset if not same as prev
            if (g_pvDBIndex.at(i).lIndexProteinFilePosition != g_pvDBIndex.at(i-1).lIndexProteinFilePosition)
            {
               g_pvProteinsList.at(g_pvProteinsList.size()-1).push_back(g_pvDBIndex.at(i).lIndexProteinFilePosition);
               g_pvDBIndex.at(i).lIndexProteinFilePosition = g_pvDBIndex.at(i-1).lIndexProteinFilePosition;
            }
         }
         else
         {
            // store new protein offset
            vector<comet_fileoffset_t> temp;
            temp.push_back(g_pvDBIndex.at(i).lIndexProteinFilePosition);
            g_pvProteinsList.push_back(temp);
            g_pvDBIndex.at(i).lIndexProteinFilePosition = lProtCount;
            lProtCount++;
         }
      }
   }

   g_pvDBIndex.erase(unique(g_pvDBIndex.begin(), g_pvDBIndex.end()), g_pvDBIndex.end());

   // sort by mass;
   sort(g_pvDBIndex.begin(), g_pvDBIndex.end(), CompareByMass);

/*
   for (std::vector<DBIndex>::iterator it = g_pvDBIndex.begin(); it != g_pvDBIndex.end(); ++it)
   {
      printf("OK after unique ");
      if ((*it).pcVarModSites[strlen((*it).szPeptide)] != 0)
         printf("n*");
      for (unsigned int x = 0; x < strlen((*it).szPeptide); x++)
      {
         printf("%c", (*it).szPeptide[x]);
         if ((*it).pcVarModSites[x] != 0)
            printf("*");
      }
      if ((*it).pcVarModSites[strlen((*it).szPeptide) + 1] != 0)
         printf("c*");
      printf("   %f   %lld\n", (*it).dPepMass, (*it).lIndexProteinFilePosition);
   }
   printf("\n");
*/

   sprintf(szOut, " - writing file\n");
   logout(szOut);
   fflush(stdout);

   // write out index header
   fprintf(fptr, "Comet indexed database.\n");
   fprintf(fptr, "InputDB:  %s\n", g_staticParams.databaseInfo.szDatabase);
   fprintf(fptr, "MassRange: %f %f\n", g_staticParams.options.dPeptideMassLow, g_staticParams.options.dPeptideMassHigh);
   fprintf(fptr, "MassType: %d %d\n", g_staticParams.massUtility.bMonoMassesParent, g_staticParams.massUtility.bMonoMassesFragment);
   fprintf(fptr, "Enzyme: %s [%d %s %s]\n", g_staticParams.enzymeInformation.szSearchEnzymeName,
      g_staticParams.enzymeInformation.iSearchEnzymeOffSet,
      g_staticParams.enzymeInformation.szSearchEnzymeBreakAA,
      g_staticParams.enzymeInformation.szSearchEnzymeNoBreakAA);
   fprintf(fptr, "NumPeptides: %ld\n", (long)g_pvDBIndex.size());

   // write out static mod params A to Z is ascii 65 to 90 then terminal mods
   fprintf(fptr, "StaticMod:");
   for (int x = 65; x <= 90; x++)
      fprintf(fptr, " %f", g_staticParams.staticModifications.pdStaticMods[x]);
   fprintf(fptr, " %f", g_staticParams.staticModifications.dAddNterminusPeptide);
   fprintf(fptr, " %f", g_staticParams.staticModifications.dAddCterminusPeptide);
   fprintf(fptr, " %f", g_staticParams.staticModifications.dAddNterminusProtein);
   fprintf(fptr, " %f\n", g_staticParams.staticModifications.dAddCterminusProtein);

   // write out variable mod params
   fprintf(fptr, "VariableMod:");
   for (int x = 0; x < VMODS; x++)
      fprintf(fptr, " %s %f", g_staticParams.variableModParameters.varModList[x].szVarModChar, g_staticParams.variableModParameters.varModList[x].dVarModMass);
   fprintf(fptr, "\n\n");

   int iTmp = (int)g_pvProteinNames.size();
   comet_fileoffset_t *lProteinIndex = new comet_fileoffset_t[iTmp];
   for (int i = 0; i < iTmp; i++)
      lProteinIndex[i] = -1;

   // first just write out protein names. Track file position of each protein name
   int ctProteinNames = 0;
   for (auto it = g_pvProteinNames.begin(); it != g_pvProteinNames.end(); ++it)
   {
      lProteinIndex[ctProteinNames] = comet_ftell(fptr);
      fwrite(it->second.szProt, sizeof(char)*WIDTH_REFERENCE, 1, fptr);
      it->second.iWhichProtein = ctProteinNames;
      ctProteinNames++;
   }

   // next write out the peptides and track peptide mass index
   int iMaxPeptideMass = (int)(g_staticParams.options.dPeptideMassHigh);
   comet_fileoffset_t *lIndex = new comet_fileoffset_t[iMaxPeptideMass];
   for (int x = 0; x < iMaxPeptideMass; x++)
      lIndex[x] = -1;

   // write out struct data
   int iPrevMass = 0;
   int iWhichProtein = 0;
   long lNumMatchedProteins = 0;
   for (std::vector<DBIndex>::iterator it = g_pvDBIndex.begin(); it != g_pvDBIndex.end(); ++it)
   {
      if ((int)((*it).dPepMass) > iPrevMass)
      {
         iPrevMass = (int)((*it).dPepMass);
         if (iPrevMass < iMaxPeptideMass)
            lIndex[iPrevMass] = comet_ftell(fptr);
      }

      fwrite(&(*it), sizeof(struct DBIndex), 1, fptr);

      // now write out all duplicate proteins file positions
      long lSize = (long)g_pvProteinsList.at((*it).lIndexProteinFilePosition).size();
      lNumMatchedProteins = lSize;
      fwrite(&lNumMatchedProteins, sizeof(long), 1, fptr);
      for (long x = 0; x < lSize; x++)
      {
         // find protein by matching g_pvProteinNames.lProteinFilePosition to g_pvProteinNames.lProteinIndex;
         auto result = g_pvProteinNames.find((g_pvProteinsList.at((*it).lIndexProteinFilePosition)).at(x));
         if (result != g_pvProteinNames.end())
         {
            iWhichProtein = result->second.iWhichProtein;
         }

         fwrite(&(lProteinIndex[iWhichProtein]), sizeof(comet_fileoffset_t), 1, fptr);
      }

   }

   comet_fileoffset_t lEndOfPeptides = comet_ftell(fptr);

   int iTmpCh = (int)(g_staticParams.options.dPeptideMassLow);
   fwrite(&iTmpCh, sizeof(int), 1, fptr);  // write min mass
   fwrite(&iMaxPeptideMass, sizeof(int), 1, fptr);  // write max mass
   size_t tNumPeptides = g_pvDBIndex.size();
   fwrite(&tNumPeptides, sizeof(size_t), 1, fptr);  // write # of peptides
   fwrite(lIndex, sizeof(comet_fileoffset_t), iMaxPeptideMass, fptr); // write index
   fwrite(&lEndOfPeptides, sizeof(comet_fileoffset_t), 1, fptr);  // write ftell position of min/max mass, # peptides, peptide index

   fclose(fptr);

   sprintf(szOut, " - done\n");
   logout(szOut);
   fflush(stdout);

   CometSearch::DeallocateMemory(g_staticParams.options.iNumThreads);

   g_pvDBIndex.clear();
   g_pvProteinNames.clear();
   delete[] lProteinIndex;
   delete[] lIndex;

   return bSucceeded;
}
