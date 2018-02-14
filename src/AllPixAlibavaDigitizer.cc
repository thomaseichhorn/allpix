/**
 *  Author:
 *    Thomas Eichhorn <thomas.eichhorn@desy.de>
 *
 *  allpix Authors:
 *   John Idarraga <idarraga@cern.ch>
 *   Mathieu Benoit <benoit@lal.in2p3.fr>
 */

#include "AllPixAlibavaDigitizer.hh"
#include "AllPixTrackerHit.hh"
#include "G4DigiManager.hh"
#include "G4VDigitizerModule.hh"
#include "AllPixGeoDsc.hh"

#include "TMath.h"
#include "TF1.h"

#include "CLHEP/Random/RandGauss.h"

#include "TFile.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TProfile.h"


AllPixAlibavaDigitizer::AllPixAlibavaDigitizer ( G4String modName, G4String hitsColName, G4String digitColName ) : AllPixDigitizerInterface ( modName )
{
    // Registration of digits collection name
    collectionName.push_back ( digitColName );
    m_hitsColName.push_back ( hitsColName );

    // threshold
    m_digitIn.thl = 0.0;

    InitVariables ( );

}


double AllPixAlibavaDigitizer::funcSharemean ( double x, double shareval, double pitch )
{
    double half_pitch = pitch / 2.0;

    double shared_percent = 0.0;
    // constant charge sharing
    //shared_percent = shareval;

    // linear charge sharing
    //shared_percent = fabs ( ( shareval / half_pitch ) * x );

    // quadratic charge sharing
    shared_percent = ( shareval / ( half_pitch * half_pitch ) ) * x * x;
    return shared_percent;
}


void AllPixAlibavaDigitizer::InitVariables ( )
{
    // set to true for some verbose output
    debugmode = false;

    AllPixGeoDsc * gD = GetDetectorGeoDscPtr ( );

    // % / 100 mean of charge to be shared, passed from the resistivity, default: 0.225
    sharemean = gD -> GetResistivity ( );

    // max charge shared at edge of cell, passed from GetSaturationEnergy
    edgemaxsharing = gD -> GetSaturationEnergy ( );

    // % / 100 sigma of charge to be shared, passed from clockunit, default: 0.05
    sharesig = gD -> GetClockUnit ( );

    // 1- x % / 100 by how much signals are scaled down if they are not in the centre of the pixel
    // 0.14132 / 44.479 1 / um from fit in run000022 epi 100p unirr data
    // 0.02712 / 28.0768 1 / um from fit in run000225 epi 100p 1.3e16 data
    // passed from mipcharge
    distancereduce = gD -> GetMIPCharge ( );

    // gain factor for signals going from electrons to alibava adcs, default 1 / 183.0
    gain = gD -> GetChipNoise ( );

    // threshold for digitizer output, passed from chipthreshold, default 0.0
    threshold = gD -> GetThreshold ( );

    // vars for output control plots, not relevant for the simulation

    // seedcut for plotting
    seedcut = 5.0;

    // clustercut for plotting
    clustercut = 2.5;

    // expected noise mean
    noisemean = 0.0;

    // expected noise sigma
    noisesigma = 5.0;

    NPixelX = gD -> GetNPixelsX ( );
    NPixelY = gD -> GetNPixelsY ( );
    pitchX = gD -> GetPixelX ( ) / um;
    pitchY = gD -> GetPixelY ( ) / um;

    // give some output
    G4cout << " " << endl;
    G4cout << "Alibava Digitizer Settings:" << endl;
    G4cout << "Charge sharing mean is " << sharemean << "!" << G4endl;
    G4cout << "Charge sharing max at edge is " << edgemaxsharing << "!" << G4endl;
    G4cout << "Charge sharing sigma is " << sharesig << "!" << G4endl;
    G4cout << "Charge deposition distance reduction is " << distancereduce << "!" << G4endl;
    G4cout << "Gain is " << gain << "!" << G4endl;
    G4cout << "Digitizer threshold is " << threshold << "!" << G4endl;

    rootoutputfile = new TFile ( "someoutputfile.root", "RECREATE" );

    histochargeleft = new TH1D ( "Output Charge Left Channel", "Output Charge Left Channel", 211, -10, 200 );
    histochargeleft -> SetTitle ( "Output Charge Left Channel;ADCs;Entries" );

    histochargemiddle = new TH1D ( "Output Charge Middle Channel", "Output Charge Middle Channel", 211, -10, 200 );
    histochargemiddle -> SetTitle ( "Output Charge Middle Channel;ADCs;Entries" );

    histochargeright = new TH1D ( "Output Charge Right Channel", "Output Charge Right Channel", 211, -10, 200 );
    histochargeright -> SetTitle ( "Output Charge Right Channel;ADCs;Entries" );

    histohitcount = new TH1I ( "Hit Count per Track", "Hit Count per Track", 11, -0.5, 10.5 );
    histohitcount -> SetTitle ( "Hit Count per Track;Hits per Track;Entries" );

    histohitcountoverthresh = new TH1I ( "Hit Count per Track over Digitizer Threshold", "Hit Count per Track over Digitizer Threshold", 11, -0.5, 10.5 );
    histohitcountoverthresh -> SetTitle ( "Hit Count per Track over Digitizer Threshold;Hits per Track;Entries" );

    histoclustersize = new TH1I ( "CShareclustersize", "CShareclustersize", 5, -0.5, 4.5 );
    histoclustersize -> SetTitle ( "Expected Cluster Size from Charge Sharing;Cluster Size;Entries" );

    histoloss_position = new TH1D ( "Charge Lost due to Track Position", "Charge Lost due to Track Position", 250, 0, 250 );
    histoloss_position -> SetTitle ( "Charge Lost due to Track Position;ADCs;Entries" );

    histoloss_share = new TH1D ( "Charge Lost due to Charge Sharing", "Charge Lost due to Charge Sharing", 250, 0, 250 );
    histoloss_share -> SetTitle ( "Charge Lost due to Charge Sharing;ADCs;Entries" );

    histoadcout = new TH1D ( "Track Induced Signal in ADCs", "Track Induced Signal in ADCs", 250, 0, 250 );
    histoadcout -> SetTitle ( "Track Induced Signal in ADCs;ADCs;Entries" );

    histoeleout = new TH1D ( "Track Induced Signal in Electrons", "Track Induced Signal in Electrons", 250, 0, 40000 );
    histoeleout -> SetTitle ( "Track Induced Signal in Electrons;Electrons;Entries" );

    histoshare_position = new TProfile ( "Charge Shared vs. Track Position", "Charge Shared vs. Track Position", 100, 0, 1 );
    histoshare_position -> SetTitle ( "Charge Shared vs. Track Position;Track Position in Units of Pitch;ADCs" );

    histohitmap_pix = new TH1I ( "Track Point Hitmap within Channel", "Track Point Hitmap within Channel", 100, 0, 1 );
    histohitmap_pix -> SetTitle ( "Track Point Hitmap within Channel;Track Position in Units of Pitch;Entries" );

    // strip sensor: assume sensitive dimesion has more pixels with smaller pitch than the other
    if ( NPixelX > NPixelY && pitchX < pitchY )
    {
	G4cout << "Assuming sensitive dimension is x!" << G4endl;
	sensdirection = "x";

	histohitmap = new TH1I ( "Track Point Hitmap", "Track Point Hitmap", NPixelX, 0, NPixelX - 1 );
	histohitmap -> SetTitle ( "Track Point Hitmap;Channel in X;Entries" );

    }
    else if ( NPixelY > NPixelX && pitchY < pitchX )
    {
	G4cout << "Assuming sensitive dimension is y!" << G4endl;
	sensdirection = "y";

	histohitmap = new TH1I ( "Track Point Hitmap", "Track Point Hitmap", NPixelY, 0, NPixelY - 1 );
	histohitmap -> SetTitle ( "Track Point Hitmap;Channel in Y;Entries" );
    }
    else
    {
	G4cout << "Could not determine sensitive dimension!" << G4endl;
	exit ( -1 );
    }

}


AllPixAlibavaDigitizer::~AllPixAlibavaDigitizer ( )
{
    rootoutputfile -> cd ( );

    TF1 * fitleft = new TF1 ( "fitleft", "landau" );
    histochargeleft -> Fit ( fitleft, "Q" );
    histochargeleft -> Write ( );

    TF1 * fitmiddle = new TF1 ( "fitmiddle", "landau" );
    histochargemiddle -> Fit ( fitmiddle, "Q" );
    histochargemiddle -> Write ( );

    TF1 * fitright = new TF1 ( "fitright", "landau" );
    histochargeright -> Fit ( fitright, "Q" );
    histochargeright -> Write ( );

    histohitcount -> Write ( );
    histohitcountoverthresh -> Write ( );
    histoclustersize -> Write ( );

    histoloss_position -> Write ( );
    histoloss_share -> Write ( );

    TF1 * fitadc = new TF1 ( "fitadc", "landau" );
    histoadcout -> Fit ( fitadc, "Q" );
    histoadcout -> Write ( );
    TF1 * fitelec = new TF1 ( "fitelec", "landau" );
    histoeleout -> Fit ( fitelec, "Q" );
    histoeleout -> Write ( );

    histohitmap -> Write ( );

    histoshare_position -> Write ( );

    histohitmap_pix -> Write ( );

    rootoutputfile -> Close ( );

}


void AllPixAlibavaDigitizer::Digitize ( )
{
    // create the digits collection
    m_digitsCollection = new AllPixAlibavaDigitsCollection ( "AllPixAlibavaDigitizer", collectionName[0] );

    // get the digiManager
    G4DigiManager * digiMan = G4DigiManager::GetDMpointer ( );

    // BoxSD_0_HitsCollection
    G4int hcID = digiMan -> GetHitsCollectionID ( m_hitsColName[0] );

    AllPixTrackerHitsCollection * hitsCollection = 0;
    hitsCollection = ( AllPixTrackerHitsCollection* ) ( digiMan -> GetHitsCollection ( hcID ) );

    // temporary data structure
    map < pair < G4int, G4int >, G4double > pixelsContent;

    // the pixel pointed to
    pair < G4int, G4int > tempPixel;

    // the left neighbour
    pair < G4int, G4int > leftcsharePixel;

    // the right neighbour
    pair < G4int, G4int > rightcsharePixel;

    // the entry count
    G4int nEntries = hitsCollection -> entries ( );

    // summed signals in the 3 channels
    G4double totalcentresignal = 0.0;
    G4double totalleftsignal = 0.0;
    G4double totalrightsignal = 0.0;

    // summed count of the charge lost due to position
    G4double totallosssignal = 0.0;

    // summed count of the charge lost due to charge sharing
    G4double totalshare = 0.0;

    // the output signal in electrons
    G4double electronsignal = 0.0;

    // the output signal in adcs
    G4double adcsignal = 0.0;

    // electrons created per eV
    G4double elec = 3.64;

    // the loop
    for ( G4int itr  = 0; itr < nEntries; itr++ )
    {

	if ( debugmode == true )
	{
	    G4cout << "Looping entries..." << endl;
	}

	// positions
	tempPixel.first  = ( *hitsCollection ) [itr] -> GetPixelNbX ( );
	tempPixel.second = ( *hitsCollection ) [itr] -> GetPixelNbY ( );

	G4ThreeVector vec = ( *hitsCollection )[itr] -> GetPosWithRespectToPixel ( ) / um;

	// the signal for the individual pixels
	G4double leftsignal[2] = { 0.0 };
	G4double centresignal[2] = { 0.0 };
	G4double rightsignal[2] = { 0.0 };

	// the signal in electrons
	G4int integerelectron = ( ( ( *hitsCollection ) [itr] -> GetEdep ( ) ) / eV ) / elec;
	electronsignal += integerelectron;

	// the signal with gain (in MeV)
	//G4double signal = ( gain* ( ( *hitsCollection ) [itr] -> GetEdep ( ) ) );
	G4double signal = gain * integerelectron;
	adcsignal += signal;

	if ( debugmode == true )
	{
	    G4cout << " " << G4endl;
	    G4cout << "Signal with gain is " << signal << " ADCs" << G4endl;
	    G4cout << "Actual deposition " << ( ( *hitsCollection ) [itr] -> GetEdep ( ) ) / MeV << " MeV = " << ( ( ( *hitsCollection ) [itr] -> GetEdep ( ) ) / eV ) / elec << " e-" << G4endl;

	    // where are we hitting the pixel?
	    G4cout << "Position wrt pixel centre is " << vec[0] << " um in x and " << vec[1] << " um in y" << G4endl;
	}

	// reduce charge depending on position relative to the pixel centre
	double positionloss[2] = { 0.0 };
	if ( sensdirection == "x" )
	{
	    positionloss[0] = distancereduce * fabs ( vec[0] ) * signal;
	    positionloss[1] = 0.0;
	}
	if ( sensdirection == "y" )
	{
	    positionloss[0] = 0.0;
	    positionloss[1] = distancereduce * fabs ( vec[1] ) * signal;
	}

	if ( debugmode == true )
	{
	    G4cout << "Position loss is " << positionloss[0] << " ADCs in x and " << positionloss[1] << " ADCs in y" << G4endl;
	}

	// the charge to deposit on the centre and to be shared to both neighbours
	double chargetoshare[2] = { 0.0 };
	if ( sensdirection == "x" )
	{
	    chargetoshare[0] = signal - positionloss[0];
	    chargetoshare[1] = 0.0;
	}
	if ( sensdirection == "y" )
	{
	    chargetoshare[0] = 0.0;
	    chargetoshare[1] = signal - positionloss[1];
	}

	if ( debugmode == true )
	{
	    G4cout << "Charge for deposition and sharing is " << chargetoshare[0] << " ADCs in x and " << chargetoshare[1] << " ADCs in y" << G4endl;
	}

	// how much is 'lost' due to the position in the pixel
	totallosssignal += positionloss[0];
	totallosssignal += positionloss[1];

	// chargeshare spreads a percentage of the actual charge to the two neighbouring pixels
	// this is randomized
	double chargeshare = 0.0;
	
	if ( sensdirection == "x" )
	{
	    chargeshare = CLHEP::RandGauss::shoot ( funcSharemean ( vec[0], edgemaxsharing, pitchX ), sharesig );
	}
	else if ( sensdirection == "y" )
	{
	    chargeshare = CLHEP::RandGauss::shoot ( funcSharemean ( vec[1], edgemaxsharing, pitchY ), sharesig );
	}

	if ( debugmode == true )
	{
	    G4cout << "Sharing " << chargeshare * 100.0 << " % of the deposition charge" << G4endl;
	}

	// for the centre
	centresignal[0] = chargetoshare[0] - ( chargetoshare[0] * chargeshare );
	centresignal[1] = chargetoshare[1] - ( chargetoshare[1] * chargeshare );
	totalshare += chargetoshare[0] * chargeshare;
	totalshare += chargetoshare[1] * chargeshare;

	// assume this is right...
	if ( vec[0] > 0.0 )
	{
	    rightsignal[0] = ( chargetoshare[0] * chargeshare ) / 2.0;
	    if ( debugmode == true )
	    {
		G4cout << "Sharing right in x: " << rightsignal[0] << G4endl;
	    }
	}
	if ( vec[0] < 0.0 )
	{
	    leftsignal[0] = ( chargetoshare[0] * chargeshare ) / 2.0;
	    if ( debugmode == true )
	    {
		G4cout << "Sharing left in x: " << leftsignal[0] << G4endl;
	    }
	}
	if ( vec[1] > 0.0 )
	{
	    rightsignal[1] = ( chargetoshare[1] * chargeshare ) / 2.0;
	    if ( debugmode == true )
	    {
		G4cout << "Sharing right in y: " << rightsignal[1] << G4endl;
	    }
	}
	if ( vec[1] < 0.0 )
	{
	    leftsignal[1] = ( chargetoshare[1] * chargeshare ) / 2.0;
	    if ( debugmode == true )
	    {
		G4cout << "Sharing left in y: " << leftsignal[1] << G4endl;
	    }
	}

	// positions left and right, catch that there is only one pixel in the unsensitive dimension
	pixelsContent[tempPixel] += ( centresignal[0] + centresignal[1] );
	if ( sensdirection == "y" )
	{
	    // catch sensor edges
	    if ( tempPixel.second > 0 && tempPixel.second < ( NPixelY - 1 ) )
	    {
		leftcsharePixel.first = tempPixel.first;
		leftcsharePixel.second = tempPixel.second - 1;
		rightcsharePixel.first = tempPixel.first;
		rightcsharePixel.second = tempPixel.second + 1;
	    }
	    else
	    {
		leftcsharePixel.first = tempPixel.first;
		leftcsharePixel.second = tempPixel.second;
		rightcsharePixel.first = tempPixel.first;
		rightcsharePixel.second = tempPixel.second;
	    }
	}
	else if ( sensdirection == "x" )
	{
	    // catch sensor edges
	    if ( tempPixel.first > 0 && tempPixel.first < ( NPixelX - 1 ) )
	    {
		leftcsharePixel.first = tempPixel.first - 1;
		leftcsharePixel.second = tempPixel.second;
		rightcsharePixel.first = tempPixel.first + 1;
		rightcsharePixel.second = tempPixel.second;
	    }
	    else
	    {
		leftcsharePixel.first = tempPixel.first;
		leftcsharePixel.second = tempPixel.second;
		rightcsharePixel.first = tempPixel.first;
		rightcsharePixel.second = tempPixel.second;
	    }
	}
	else
	{
	    leftcsharePixel.first = tempPixel.first - 1;
	    leftcsharePixel.second = tempPixel.second - 1;
	    rightcsharePixel.first = tempPixel.first + 1;
	    rightcsharePixel.second = tempPixel.second + 1;
	}
	pixelsContent[leftcsharePixel] += ( leftsignal[0] + leftsignal[1] );
	pixelsContent[rightcsharePixel] += ( rightsignal[0] + rightsignal[1] );

	totalleftsignal += ( leftsignal[0] + leftsignal[1] );
	totalcentresignal += ( centresignal[0] + centresignal[1] );
	totalrightsignal += ( rightsignal[0] + rightsignal[1] );

	if ( debugmode == true )
	{
	    G4cout << "Left   pixel is ( " << leftcsharePixel.first << " | " << leftcsharePixel.second << " ) with charge " << pixelsContent[leftcsharePixel] << G4endl;
	    G4cout << "Middle pixel is ( " << tempPixel.first << " | " << tempPixel.second << " ) with charge " << pixelsContent[tempPixel] << G4endl;
	    G4cout << "Right  pixel is ( " << rightcsharePixel.first << " | " << rightcsharePixel.second << " ) with charge " << pixelsContent[rightcsharePixel] << G4endl;
	}

	if ( sensdirection == "x" )
	{
	    double relpos = vec[0] / pitchX;
	    if ( relpos < 0 )
	    {
		relpos += 1.0;
	    }
	    histoshare_position -> Fill ( relpos, ( rightsignal[0] + leftsignal[0] ) );
	    histohitmap_pix -> Fill ( relpos );
	}
	else if ( sensdirection == "y" )
	{
	    double relpos = vec[1] / pitchY;
	    if ( relpos < 0 )
	    {
		relpos += 1.0;
	    }
	    histoshare_position -> Fill ( relpos, ( rightsignal[1] + leftsignal[1] ) );
	    histohitmap_pix -> Fill ( relpos );
	}

    }

    if ( debugmode == true )
    {
	G4cout << " " << G4endl;
	G4cout << "Total Track output signal is " << adcsignal << " ADCs, or " << electronsignal << " e-" << G4endl;
    }

    if ( adcsignal > 0.0 )
    {
	histoadcout -> Fill ( adcsignal );
    }
    if ( electronsignal > 0.0 )
    {
	histoeleout -> Fill ( electronsignal );
    }
    if ( totalleftsignal > 0.0 )
    {
	histochargeleft -> Fill ( totalleftsignal );
    }
    if ( totalcentresignal > 0.0 )
    {
	histochargemiddle -> Fill ( totalcentresignal );
    }
    if ( totalrightsignal > 0.0 )
    {
	histochargeright -> Fill ( totalrightsignal );
    }
    histoloss_position -> Fill ( totallosssignal );
    histoloss_share -> Fill ( totalshare );



    double noise = CLHEP::RandGauss::shoot ( noisemean, noisesigma );

    int clustersize = 0;
    if ( ( totalcentresignal + noise ) > seedcut * noisesigma )
    {
	clustersize++;

	noise = CLHEP::RandGauss::shoot ( noisemean, noisesigma );
	if ( ( totalleftsignal + noise ) > clustercut * noisesigma )
	{
	    clustersize++;
	}
	noise = CLHEP::RandGauss::shoot ( noisemean, noisesigma );
	if ( ( totalrightsignal + noise ) > clustercut * noisesigma )
	{
	    clustersize++;
	}
    }
    histoclustersize -> Fill ( clustersize );

    // Now create digits, one per pixel
    // Second entry in the map is the energy deposit in the pixel
    map < pair < G4int, G4int >, G4double >::iterator pCItr = pixelsContent.begin ( );

    // NOTE that there is a nice interface which provides useful info for hits.
    // For instance, the following member gives you the position of a hit with respect
    //  to the center of the pixel.
    // G4ThreeVector vec = (*hitsCollection)[itr]->GetPosWithRespectToPixel();
    // See class AllPixTrackerHit !

    // Also, you have full access to the detector geometry from this scope
    // AllPixGeoDsc * GetDetectorGeoDscPtr()
    // provides you with a pointer to the geometry description.
    // See class AllPixGeoDsc !
    G4int hitcount = 0;
    G4int hitcountoverthresh = 0;

    for ( ; pCItr != pixelsContent.end ( ); pCItr++ )
    {
	hitcount++;
	if ( ( *pCItr ).second > threshold ) // over threshold !
	{
	    hitcountoverthresh++;
	    G4double signal = 0.0;
	    AllPixAlibavaDigit * digit = new AllPixAlibavaDigit;
	    digit -> SetPixelIDX ( ( *pCItr ).first.first );
	    digit -> SetPixelIDY ( ( *pCItr ).first.second );
	    if ( sensdirection == "x" )
	    {
		histohitmap -> Fill ( ( *pCItr ).first.first, 1 );
	    }
	    else if ( sensdirection == "y" )
	    {
		histohitmap -> Fill ( ( *pCItr ).first.second, 1 );
	    }
	    signal = ( *pCItr ).second;
	    digit -> SetPixelCounts ( signal );
	    m_digitsCollection -> insert ( digit );
	}
    }

    histohitcount -> Fill ( hitcount );
    histohitcountoverthresh -> Fill ( hitcountoverthresh );

    G4int dc_entries = m_digitsCollection -> entries ( );
    if ( dc_entries > 0 )
    {
	G4cout << "--------> Digits Collection : " << collectionName[0] << "(" << m_hitsColName[0] << ")" << " contains " << dc_entries << " digits" << G4endl;
    }

    StoreDigiCollection ( m_digitsCollection );

}
