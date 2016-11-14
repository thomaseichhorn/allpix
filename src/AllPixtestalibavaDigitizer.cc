/**
 *  Author:
 *    thomas eichhorn <thomas.eichhorn@desy.de>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#include "AllPixPrimaryGeneratorMessenger.hh"
#include "AllPixDetectorMessenger.hh"
#include "AllPixtestalibavaDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"

#include "CLHEP/Random/RandGauss.h"

AllPixtestalibavaDigitizer::AllPixtestalibavaDigitizer(G4String modName, G4String hitsColName, G4String digitColName) 
: AllPixDigitizerInterface (modName) {

	// Registration of digits collection name
	collectionName.push_back(digitColName);
	m_hitsColName.push_back(hitsColName);

	// threshold
	m_digitIn.thl = 0.;

	InitVariables();
}

void AllPixtestalibavaDigitizer::InitVariables()
{

	double frames = 2.0;

	noisemean = 0.0;
	noisesig = 5.0/frames;
	pedestalmean = 500.0/frames;
	pedestalsig = 10.0/frames;
	commonmodemean = 0.0;
	commonmodesig = 20.0/frames;

	sharemean = 0.25;
	sharesig = 0.1;
	distancereduce = 2.0;

	AllPixGeoDsc * gD = GetDetectorGeoDscPtr();
	NPixelX = gD->GetNPixelsX();
	NPixelY = gD->GetNPixelsY();
	pitchX=gD->GetPixelX()/um;
	pitchY=gD->GetPixelY()/um;
}

AllPixtestalibavaDigitizer::~AllPixtestalibavaDigitizer(){

}

void AllPixtestalibavaDigitizer::Digitize(){

	// create the digits collection
	m_digitsCollection = new AllPixtestalibavaDigitsCollection("AllPixtestalibavaDigitizer", collectionName[0] );

	// get the digiManager
	G4DigiManager * digiMan = G4DigiManager::GetDMpointer();

	// BoxSD_0_HitsCollection
	G4int hcID = digiMan->GetHitsCollectionID(m_hitsColName[0]);

	AllPixTrackerHitsCollection * hitsCollection = 0;
	hitsCollection = (AllPixTrackerHitsCollection*)(digiMan->GetHitsCollection(hcID));

	// temporary data structure
	map<pair<G4int, G4int>, G4double > pixelsContent;
	pair<G4int, G4int> tempPixel;
	pair<G4int, G4int> leftcsharePixel;
	pair<G4int, G4int> rightcsharePixel;

	G4int nEntries = hitsCollection->entries();

	for(G4int itr  = 0 ; itr < nEntries ; itr++) {

		tempPixel.first  = (*hitsCollection)[itr]->GetPixelNbX();
		tempPixel.second = (*hitsCollection)[itr]->GetPixelNbY();
		G4ThreeVector vec = (*hitsCollection)[itr]->GetPosWithRespectToPixel()/um;

		// only y for now
		// scalefactor to reduce charge depending on position relative to the pixel centre
		// chargeshare spreads charge to the two neighbouring pixels
		double scalefactor = (pitchY-distancereduce*fabs(vec[1]))/pitchY;
		double chargeshare = 1.0 - CLHEP::RandGauss::shoot(sharemean, sharesig);
		//cout << "chargeshare " << chargeshare << endl;
		double signal = ((*hitsCollection)[itr]->GetEdep())/keV;
		signal = signal * scalefactor;
		double sharedsignal = signal * chargeshare;
		double distribsig = (signal - sharedsignal)/2.0;

		pixelsContent[tempPixel] += sharedsignal;
		leftcsharePixel.first = tempPixel.first;
		leftcsharePixel.second = tempPixel.second - 1;
		rightcsharePixel.first = tempPixel.first;
		rightcsharePixel.second = tempPixel.second + 1;
		pixelsContent[leftcsharePixel] += distribsig;
		pixelsContent[rightcsharePixel] += distribsig;
		//cout << "left pix is " << leftcsharePixel.second << " right " << rightcsharePixel.second << " sig is "<< distribsig << endl;
	}

	// Now create digits, one per pixel
	// Second entry in the map is the energy deposit in the pixel
	map<pair<G4int, G4int>, G4double >::iterator pCItr = pixelsContent.begin();

	// NOTE that there is a nice interface which provides useful info for hits.
	// For instance, the following member gives you the position of a hit with respect
	//  to the center of the pixel.
	// G4ThreeVector vec = (*hitsCollection)[itr]->GetPosWithRespectToPixel();
	// See class AllPixTrackerHit !

	// Also, you have full access to the detector geometry from this scope
	// AllPixGeoDsc * GetDetectorGeoDscPtr()
	// provides you with a pointer to the geometry description.
	// See class AllPixGeoDsc !

	// x here the unsensitive coordinate
	for (Size_t i = 0; i<NPixelX; i++)
	{
		// common mode same for all channels in this event
		G4int commonmode = CLHEP::RandGauss::shoot(commonmodemean, commonmodesig);
		for (Size_t j = 0; j<NPixelY; j++)
		{
			//cout << "pix " << i << " " << j << endl;
			AllPixtestalibavaDigit * digit = new AllPixtestalibavaDigit;
			digit->SetPixelIDX(i);
			digit->SetPixelIDY(j);
			G4int signal = 0.0;

			// add pedestal
			//signal=CLHEP::RandGauss::shoot(pedestalmean, pedestalsig);
			//cout << "pedestal " << signal << endl;

			// add noise
			//signal += CLHEP::RandGauss::shoot(noisemean, noisesig);

			// loop the channels with a hit
			for(pCItr = pixelsContent.begin() ; pCItr != pixelsContent.end() ; pCItr++)
			{

				if (((*pCItr).first.first) == i)
				{
					if (((*pCItr).first.second) == j)
					{
						if((*pCItr).second > 0) // over threshold !
						{
							signal += (*pCItr).second;
							//cout << "Particle singal on channel " << j << " is " << signal << endl;

						}
					} else {
					     //cout << " ((*pCItr).first.second) " << ((*pCItr).first.second) << " j is " << j << endl;
					}
				} else {
				    //cout << " ((*pCItr).first.first) " << ((*pCItr).first.first) << " i is " << i << endl;
				}

			}
			//signal += commonmode;
			//cout << "chan "<< j << "signal " << signal << endl;
			digit->SetPixelCounts(signal);
			m_digitsCollection->insert(digit);
		}
	}

	G4int dc_entries = m_digitsCollection->entries();
	if(dc_entries > 0){
		G4cout << "--------> Digits Collection : " << collectionName[0]
		                                           << "(" << m_hitsColName[0] << ")"
		                                           << " contains " << dc_entries
		                                           << " digits" << G4endl;
	}

	StoreDigiCollection(m_digitsCollection);

}
