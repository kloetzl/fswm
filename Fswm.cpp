/**
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details at
 * http://www.gnu.org/copyleft/gpl.html
 *
 */

#include <iostream> 
#include <stdlib.h>
#include <getopt.h>
#include <iomanip> 
#include "Sequence.h"
#include "Seed.h"
 
void printHelp(){
 std::string help = 
    "\nUsage: ./fswm [options] <sequence> "
    "\n"
    "\n<sequence> format:"  
    "\n\t Sequence must be in FASTA format. All genomes must be contained in one FASTA file. Example:"
    "\n\t >Genome1"
    "\n\t ATAGTAGATGAT.."
    "\n\t >Genome2"
    "\n\t ATAGTAGATGAT.."
    "\n\t >Genome3"
    "\n\t ATGATGATGATGATG.."
    "\n\t .."
    "\n\t "
    "\nOptions:"        
    "\n\t -h: print this help and exit"
    "\n\t -k <integer>: pattern weight (default 12)"
    "\n\t -t <integer>: numer of threads (default: 10)"
    "\n\t -s <integer>: the minimum score of a spaced-word match to be considered homologous (default: 0)"
    "\n";
	std::cout << help << std::endl;
}

int weight = 12;
int dontCare = 100;
int threads = 10;
int threshold = 0;

static std::string output_file_name = "DMat";

void parseParameters(int argc, char *argv[]){
	int option_char;
	 while ((option_char = getopt (argc, argv, "k:t:hs:o:")) != -1){
		switch (option_char){ 
			case 's': 
				threshold = atoi (optarg); 
				break;
			case 'k': 
				weight = atoi (optarg); 
				if(weight<8 || weight > 16){
					std::cerr << "Weight (-k) must be between 8 and 16"<< std::endl;
					exit (EXIT_FAILURE);
				}
				break;
			case 't': 
				threads = atoi (optarg); 
				if(threads<1){
					std::cerr << "threads (-t) must be an integer larger than 0"<< std::endl;
					exit (EXIT_FAILURE);
				}
				break;
			case 'h': 
				printHelp();
				exit (EXIT_SUCCESS);
			case 'o':
				output_file_name = optarg;
				break;
			case '?': 
				printHelp();		
				exit (EXIT_FAILURE);
      	}
	}
}

void writeDmat(std::vector<std::vector<double> > dmat, std::vector<Sequence>& sequences){
	std::ofstream regular_outfile;
	std::ostream &outfile = output_file_name == "-" ? std::cout : regular_outfile;
	if (output_file_name != "-") {
		regular_outfile.open(output_file_name);
	}

	outfile << sequences.size() << std::endl;
	for (int i = 0; i < sequences.size(); i++) 
	{
		std::string name = sequences[i].getHeader();
		for(int k = 0; k < 10; k++){
			if(k >= name.length())
				outfile << " ";
			else
				outfile << name[k];
		}
		outfile << " ";
     	for (int j = 0; j < sequences.size(); j++) 
     	{
			if (i > j) 
	    			outfile << std::fixed <<std::setprecision(12) << dmat[i][j] << "  ";
			else if(j>i)
				outfile << std::fixed<< std::setprecision(12) << dmat[j][i] << "  ";
			else
					outfile << std::setprecision(12) << "0" << "  ";
     	}
      		outfile << std::endl;
	}
	
	//std::cout<< dmat[1][0] << std::endl;
	if (output_file_name != "-") {
		regular_outfile.close();
	}
}

int main(int argc, char *argv[]){
	if(argc < 2)
	{
		printHelp();		
		exit (EXIT_FAILURE);
	}
	parseParameters(argc,  argv);
	std::string fileName(argv[argc-1]);
	Seed seed(weight,dontCare);
	std::vector<Sequence> sequences = Sequence::read(fileName);
	std::vector<std::vector<double> >DMat(sequences.size(), std::vector<double>(sequences.size(),0));
	Seed::init();
	omp_set_dynamic(0);     
	omp_set_num_threads(threads);
	//omp_set_nested(1);
	std::cerr << sequences.size() << " sequences read"<< std::endl;
	if(sequences.size() < 2){
		std::cerr << "there must be at least 2 sequences"<< std::endl;
		exit (EXIT_FAILURE);
	}
	std::cerr << "start sorting"<< std::endl;
	#pragma omp parallel for schedule(runtime)
	for(int i = 0; i < sequences.size();i++)
	{
		if(sequences[i].getSequence().size()<1000){
			std::cerr << "each sequence must be longer than 1000 base pairs"<< std::endl;
			exit (EXIT_FAILURE);
		}
		sequences[i].sortFirstBits(seed);
		sequences[i].sortFirstBitsRev(seed);
	}

	for(int i = 0; i < sequences.size(); i++)
	{
		sequences[i].sortNextBits(seed);
		sequences[i].sortNextBitsRev(seed);
	}
	std::cerr << "starting pairwise distance calculation"<< std::endl;
	int number = (sequences.size()*(sequences.size()-1))/2;
	int cnt = 1;
	//#pragma omp parallel for  num_threads(10)
	for(int i = 0; i < sequences.size(); i++)
	{
		for(int j = i + 1; j < sequences.size(); j++)
		{
			std::cerr << cnt++ << "/" << number << "\r";
			std::cerr.flush();
			DMat[i][j] = sequences[i].compareSequences(sequences[j], seed, threads, threshold);
			DMat[j][i] = DMat[i][j];
		}
	}
	std::cerr << std::endl << "done" << std::endl;
	writeDmat(DMat, sequences);
	std::cerr << std::endl << "Distances written to file " << output_file_name << std::endl;
}
