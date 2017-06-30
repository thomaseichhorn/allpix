/**
 * Author:
 *    Thomas Eichhorn <thomas.eichhorn@desy.de>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <Mathieu.Benoit@cern.ch>
 */

#ifndef AllPixAlibavaDigitizer_h
#define AllPixAlibavaDigitizer_h 1

// allpix Interface
#include "AllPixDigitizerInterface.hh"
// digits for this digitizer
#include "AllPixAlibavaDigit.hh"
#include "G4PrimaryVertex.hh"

#include <map>
#include <vector>

#include "TFile.h"
#include "TH1D.h"
#include "TH1I.h"

using namespace std;

class AllPixAlibavaDigitizer : public  AllPixDigitizerInterface
{
    public:

	AllPixAlibavaDigitizer ( G4String, G4String, G4String );
	virtual ~AllPixAlibavaDigitizer ( );

	void SetPrimaryVertex ( G4PrimaryVertex * pv )
	{
	    m_primaryVertex = pv;
	};
	void Digitize ( );
	void SetDetectorDigitInputs ( G4double )
	{

	};

    private:

	// digitInput typedef is defined in AllPixDigitizerInterface.hh
	digitInput m_digitIn;

	void InitVariables ( );

	bool debugmode;

	G4double clustercut;
	G4double distancereduce;
	G4double gain;
	G4double noisemean;
	G4double noisesigma;
	G4double pitchX;
	G4double pitchY;
	G4double seedcut;
	G4double sharemean;
	G4double sharesig;
	G4double threshold;

	int NPixelX;
	int NPixelY;

	G4String sensdirection;

	TFile * rootoutputfile;

	TH1D * histochargeleft;
	TH1D * histochargemiddle;
	TH1D * histochargeright;
	TH1D * histoloss_share;
	TH1D * histoloss_position;
	TH1D * histoeleout;
	TH1D * histoadcout;

	TH1I * histohitcount;
	TH1I * histohitcountoverthresh;
	TH1I * histoclustersize;

	AllPixAlibavaDigitsCollection * m_digitsCollection;

	vector < G4String > m_hitsColName;

	// information from EventAction
	G4PrimaryVertex * m_primaryVertex;

};

#endif
