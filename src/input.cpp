//#include <stdio.h>
//#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "input.h"

#define BUFFSZ 50000

//Set to point at the beginning of each parameter in buff
//Substitute commas or semicolon with \0
int Input::GetPars() {

	int len = strlen(buff);
	numofpars = 0;
	pars[numofpars] = buff; //First par

	for (int i=0; i < len; i++) {

		if (buff[i] == ',') { //Middle par
			buff[i] = '\0';
			sscanf( pars[numofpars], " %lf", rpars+numofpars); //Try to read par. as double

			numofpars++;
			pars[numofpars] = buff+i+1; //Next par

			continue;
		}

		if (buff[i] == ';') { //end
			buff[i] = '\0';
			sscanf( pars[numofpars], " %lf", rpars+numofpars); //Try to read par. as double

			numofpars++;
			return numofpars;
		}
	}
    return numofpars; // JEV warning
}



Input::Input(char *datafile, int emaxnumofcurves, int maxm, std::string ccdir) {

	//Open the array to hold all c. curves
	maxnumofcurves = emaxnumofcurves;
	curves = new Cal * [maxnumofcurves];
	numofcurves = 0;

	//Open the object to hold all determinations
	dets = new Dets(maxm);
	Det *tmpdet;

	//Open the array of pointers to point at each parameter
	maxnumofpars = MAXNUMOFPARS;
	pars = new char * [maxnumofpars];
	numofpars = 0;
	rpars = new double[maxnumofpars];

	plum = 0; //Flag to see if plum is being used
	int more_pars = 0; //additional parameters

	//Open a double array of doubles to hold the hiatuses locations, if any,
	//0=h's 1=alpha's 2=betas's 3=ha's and 4=hb's for each hiatus section
	hiatus_pars = new double * [5]; //pointers to rows
	for (int i=0; i<5; i++)
		hiatus_pars[i] = new double[MAXNUMOFHIATUS]; //open each row
	H = 0; //set to zero hiatuses at the beginning

	FILE *F;

  if ((F = fopen( datafile, "r")) == NULL){
        Rprintf("Could not open %s for reading\n", datafile);
        Rcpp::stop("Could not open %s for reading\n", datafile);
//		exit(-1);
	}

  Rprintf("Reading %s\n", datafile);
	char line[BUFFSZ];
	char key[10];
	int i=0, nm, j;


	do
	{
    /*	if(fgets( line, BUFFSZ, F));
		i++;
        JEV warning
*/
        if(fgets( line, BUFFSZ, F))
        i++;

		//Remove leading blank spaces
		j = 0;
		while (line[j] == ' ')
			j++ ;


		if ((line[j] == '#') || (line[j] == '\n')) //Comment or blank line, ignore line
			continue;

		//This is the syntax, key number : parameters
		if (sscanf( line, " %s %d :", key, &nm) < 2)
		{
            Rprintf("%s:%d Syntax error\n\n", datafile, i);

			break;
		}

		buff = strstr( line, ":") + 1;
		GetPars();

		/**** Debug ***
    //	printf("%d: KEY: %s, n= %d, %d pars:", i, key, nm, numofpars);//Debug
		for (int i=0; i<numofpars; i++)
        //	printf("|%s:%g|", pars[i], rpars[i]);
    //	printf("\n");
        ***************/

		//This is the calibration curve Key, we load a calibration curve here
//		double th; JEV warning
        if (strcmp( key, "Cal") == 0)
		{
			sscanf( pars[0], " %s", line); //c. curve name

			if (strcmp( "IntCal20", line) == 0) {   //int bomb
                curves[numofcurves++] = new IntCal20((int) rpars[1], ccdir);

				continue;
			}

			if (strcmp( "Marine20", line) == 0) {
                curves[numofcurves++] = new Marine20(ccdir);

				continue;
			}

			if (strcmp( "SHCal20", line) == 0) {
                curves[numofcurves++] = new SHCal20((int) rpars[1], ccdir);

				continue;
			}

			if (strcmp( "GenericCal", line) == 0) {   //int bomb
                curves[numofcurves++] = new GenericCal(pars[1]+1, ccdir); //a space after the comma

				continue;
			}

			if (strcmp( "Plum", line) == 0) {

				//                                alPhi     mPhi      alS       mS         Al
        curves[numofcurves++] = new Plum( rpars[1], rpars[2], rpars[3], rpars[4], rpars[5],
				                                  rpars[6], (int) rpars[7], pars[8]+1, ccdir); //a space after the comma
				//                                theta0,   Radon case,     supported data file,
				plum = 1; //There are 201Pb data and Plum needs to be used, flag activated.
				more_pars = ((Plum*) curves[numofcurves-1])->MorePars();

				continue;
			}

			if (strcmp( "ConstCal", line) == 0) {
				curves[numofcurves++] = new ConstCal();
				continue;
			}

            REprintf("Bacon: ERROR: Don't know how to handle curve %s, use GenericCal.\n", line);
            Rcpp::stop("Bacon: ERROR: Don't know how to handle curve %s, use GenericCal.\n", line);

//            printf("Bacon: ERROR: Don't know how to handle curve %s, use GenericCal.\n", line, line);
    //		exit(0);
		}

		if (strcmp( key, "Det") == 0){
			sscanf( pars[0], " %s", line);
					   //Det(char *enm, double ey, double estd, double x, double edeltaR, double edeltaSTD, double ea, double eb, Cal *ecc)
			tmpdet = new Det(  line   ,  rpars[1],    rpars[2], rpars[3],       rpars[4],         rpars[5],  rpars[6],  rpars[7], curves[(int) rpars[8]]);

			dets->AddDet(tmpdet);

			continue;
		}


		if (strcmp( key, "Hiatus") == 0)
		{
			hiatus_pars[0][H] = rpars[0]; //Hiatus position

			hiatus_pars[1][H] = rpars[1]; //alpha
			hiatus_pars[2][H] = rpars[2]; //beta

			hiatus_pars[3][H] = rpars[3]; //ha
			hiatus_pars[4][H] = rpars[4]; //hb

      Rprintf("Hiatus at: %f\n", hiatus_pars[0][H]);

			H++;

			//If H==0 then hiatus_pars will be ignored in the constructors BaconFix and BaconMov
			continue;
		}

		if (strcmp( key, "Bacon") == 0)
		{
			//Include the parameters for the last (first) section, or for the whole core if H == 0
			hiatus_pars[0][H] = -10.0; //Hiatus position, in this case, this one will be ignored

			hiatus_pars[1][H] = rpars[7+1]; //alpha
			hiatus_pars[2][H] = rpars[7+2]; //beta

			hiatus_pars[3][H] = 0.0; //ha, in this case, this one will be ignored
			hiatus_pars[4][H] = 0.0; //hb, in this case, this one will be ignored

			unsigned long int seed;
			if (numofpars == 11) // should this be 12?
				seed = 0; //automatic seed set with time()
			else
				seed = (unsigned long int) rpars[12];

			//Open the Bacon object
			sscanf( pars[0], " %s", line); //type of object

			if (strcmp( line, "FixNor") == 0){
					                //dets                K   H  hiatus_pars         a         b
				bacon = new BaconFix( dets,  (int) rpars[1],  H, hiatus_pars, rpars[6], rpars[7],
					rpars[2], rpars[3], rpars[4], rpars[5], rpars[10], rpars[11], 0,      seed, more_pars);
					//MinYr     MaxYr       th0      thp0         c0       cm   useNor
			}

			if (strcmp( line, "FixT") == 0){
					                //dets                K   H  hiatus_pars         a         b
				bacon = new BaconFix( dets,  (int) rpars[1],  H, hiatus_pars, rpars[6], rpars[7],
					rpars[2], rpars[3], rpars[4], rpars[5], rpars[10], rpars[11], 1,      seed, more_pars);
					//MinYr     MaxYr       th0      thp0         c0       cm  useT
			}


			//Then open the twalk object
			BaconTwalk = new twalk( *bacon, bacon->Getx0(), bacon->Getxp0(), bacon->get_dim());

			break;  //this should be the last key in the program
		}

        Rprintf("Unknown key: %s\n", key);

	} while (!feof(F));

	bacon->ShowDescrip();

  Rprintf("\n");
}

void Input::outputFiles(std::string outputfile1){
	if(! isPlum() ) return; //the output file is correct

	//printf("Plum is needed\n");
	FILE *F, *F1, *F2/*, *F3*/;
	char line[BUFFSZ];
	const char delim[3] = "\t ";
	std::string delimiter = ".out";
	char *token;
	int lin = 0;

	double *x, w;
	int K;

	if ((F = fopen( outputfile1.c_str(), "r")) == NULL){
		Rprintf("Could not open %s for reading\n", outputfile1.c_str());
		Rcpp::stop("Could not open %s for reading\n", outputfile1.c_str());
		//exit(-1);
	}
	//int BUFFSZ = 1024;

	std::string split = outputfile1;

	std::string outputFilePlum1 = split.substr(0, split.find(delimiter)) + "_bacon.out";
	std::string outputFilePlum2 = split.substr(0, split.find(delimiter)) + "_plum.out";

	//std::string outputFilePlum3 = split.substr(0, split.find(delimiter)) + "_PyPlum.csv";


	if ((F1 = fopen( outputFilePlum1.c_str(), "w" )) == NULL){
		Rprintf("Could not open %s for writing\n", outputFilePlum1.c_str());
		Rcpp::stop("Could not open %s for writing\n", outputFilePlum1.c_str());
		//exit(-1);
	}

	if ((F2 = fopen( outputFilePlum2.c_str(), "w" )) == NULL){
		Rprintf("Could not open %s for writing\n", outputFilePlum2.c_str());
		Rcpp::stop("Could not open %s for writing\n", outputFilePlum2.c_str());
		//exit(-1);
	}

	/*if ((F3 = fopen( outputFilePlum3.c_str(), "w" )) == NULL){
		Rprintf("Could not open %s for writing\n", outputFilePlum3.c_str());
		Rcpp::stop("Could not open %s for writing\n", outputFilePlum3.c_str());
		//exit(-1);
	}*/



	while(fgets( line, BUFFSZ, F) != NULL){
		lin++;
		if(lin==1) continue; //jump the first blank line

		token = strtok(line, delim);
		std::vector<double> vd;

		while( token != NULL ) {
			double tmp = atof(token);
			vd.push_back(tmp);
			token = strtok(NULL, delim);
		}

		for (unsigned long int i = vd.size() - (GetnPs()+2) ; i < vd.size()-1; i++) { // was just int i; MB 27 Jan 2020
			fprintf(F2, "\t%13.6g", vd[i]); //Writes the rest of the lines to _plum.out
		}
		fprintf(F2, "\n");
		//fprintf(F3, "\n"); dont change the line, the rest of the lines will be added next

		K = vd.size() - ( GetnPs()+2 ) - 2;
		x  = new double[K+2];
		x[0] = vd[0]; //th0
		x[K+1] = vd[K+1]; //w
		w = x[K+1];
		x[K]  = vd[K]; //= alpha[K]
		for (int k=K-1; k>0; k--) {
			x[k]  = w*x[k+1] + (1.0-w)*vd[k]; //Create the x's
		}
		//WATCH OUT ... NO HIATUSES!!!!!!

		for (int i = 0; i < K+2; i++) {
			fprintf(F1, "\t%13.6g", x[i]); //Writes each line to _bacon.out
		}
		fprintf(F1, "\t%13.6g\n", vd[ vd.size()-1 ]); //and adds the energy

		/*fprintf(F3, " %13.6g, ", x[K+1]); //Writes the w to _PyPlum
		for (int i = 1; i < K+1; i++) {
			fprintf(F3, " %13.6g, ", vd[i]); //Writes the alphas to _PyPlum ... theta0 is ignored
		}
		fprintf(F3, " %13.6g\n", vd[ vd.size()-1 ]); //and adds the energy to _PyPlum*/


		// delete x; MB commented this March 25 2022 to get rid of gcc debian warning
		delete[] x;

	}

	fclose(F);
	fclose(F1);
	fclose(F2);

	if(remove( outputfile1.c_str()) !=0 ){
		REprintf("PLUM: ERROR: Couldn't remove the file %s\n", outputfile1.c_str());
		Rcpp::stop("PLUM: ERROR: Couldn't remove the file %s\n", outputfile1.c_str());
	}
	if(rename( outputFilePlum1.c_str(), outputfile1.c_str() ) !=0 ){
		REprintf("PLUM: ERROR: Couldn't create the file %s\n", outputfile1.c_str());
		Rcpp::stop("PLUM: ERROR: Couldn't create the file %s\n", outputfile1.c_str());
	}

}
