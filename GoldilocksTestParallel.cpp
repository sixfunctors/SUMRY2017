// GoldilocksTestParallel.cpp
// by Connor Halleck-Dube
// Adapted from code ("usefcns.h" and "functions.cpp") by Nicolle Gruzling
// Developed as part of SUMRY 2017
// https://arxiv.org/abs/1709.03663

// This is the second part of an algorithm for the enumeration of Goldilocks
// linear threshold functions (GLTFs). This program recieves a list of candidate 
// generators for the set of Goldilocks functions, and tests them for 
// separability using an array of tester threads. They then pass the number of
// Goldilocks functions in the orbit of each generator to a totaler thread, 
// which outputs the total number of Goldilocks functions on n variables.

//					     /--> tester[i-1] --\
// Main -----> candq -------> tester[i]   ---------> countq ---------> totaler
//	 	  [candidates]   \--> tester[i+1] --/   [partial counts]	


#include "usefcns.h"
#include "bigint.h"
#include <fstream>
#include <mutex>
#include <thread>
#include <condition_variable>
#include <chrono>
#include <bitset>
#include <sstream>
#include <new>
#include <vector>
#include <iostream>

// Program uses concurrentqueue, an implementation of a thread-safe
// blocking queue. Concurrentqueue is published under the Simplified 
// BSD License. 
// Copyright (c) 2013-2016, Cameron Desrochers. All rights reserved.
#include "blockingconcurrentqueue.h"

using namespace std; 


// WARNING: For ease of parallelization and error-checking, 
// this version of the program requires the total # of LTFs on n variables 
// TOTALT to be already known. 

//const unsigned n = 3; const unsigned tn = 8;
//const unsigned n=4; const unsigned tn = 16;
//const unsigned n=5; const unsigned tn = 32; const int TOTALT = 21;
//const unsigned n=6; const unsigned tn = 64; const int TOTALT = 135;
//const unsigned n=7; const unsigned tn = 128; const int TOTALT = 2470;
//const unsigned n=8; const unsigned tn = 256; const int TOTALT = 319124;
//
const unsigned n=9; const unsigned tn = 512; const int TOTALT = 1214554343;

vector<int> Great[tn];
vector<int> Less[tn];

#include "functions.cpp"

// Two thread-safe queues. One for the functions, another for the sums
moodycamel::BlockingConcurrentQueue<bitset<tn>> candq;
moodycamel::BlockingConcurrentQueue<std::tuple<lint, int, lint, int>> countq;

// Number of threads allowed (must agree with cluster allowance)
int MAXTHREADS = 16;

// Name of the file holding the candidates
char readname[] = "/home/fas/payne_sam/cjh69/project/GoldCands9.dat";

// Name of the file holding the results
char outname[] = "/home/fas/payne_sam/cjh69/GoldCounts9.txt";

// Name of the log file
char logname[] = "/home/fas/payne_sam/cjh69/GoldLog9.txt";

// Approximate maximum number of elements in the test queue at once (loose)
int QUEUEMAX = 5000;

// How long main should wait for the queue to empty (ms)
int WAITFOR = 5;

int buflen;// Read buffer length in bitsets
const int bufsize = 2097152; // Read buffer size (in chars) (must be mult of tn)
// Smaller version 65536
char buffer[bufsize];

// Allows all threads to write to log file safely
std::mutex logMut;
void log(std::string const& msg) {
	std::lock_guard<std::mutex> lock(logMut);

	ofstream outfile;
	outfile.open(logname, ios::app);
	outfile << msg;
	outfile.close();
}

// For output writing
std::mutex outMut;
void output(std::bitset<tn> F, std::tuple<lint, lint, lint, lint>* res) {
	std::lock_guard<std::mutex> lock(outMut);

	ofstream outfile;
	outfile.open(outname, ios::app);

	outfile << F.to_string() << ";";
	outfile << std::get<0>(*res) << ",";
	outfile << std::get<1>(*res) << ",";
	outfile << std::get<2>(*res) << ",";
	outfile << std::get<3>(*res) << endl;
	outfile.close();
}

void output(std::string S) {
	std::lock_guard<std::mutex> lock(outMut);

	ofstream outfile;
	outfile.open(outname, ios::app);
	outfile << S;
	outfile.close();
}

// Thread function: Tests boolean functions F for separability. 
// If separable, puts (m, n, s, t) in retvals
// where m is the number of goldilocks functions associated to F 
// 		 n is the number of goldilocks LTFs up to symmetry associated to F
// 		 s is the number of positive, small LTFs associated to F
// 		 t is the number of PS LTFs up to symmetry associated to F
// Passes results into a shared queue, where they are combined
void tester(int id){

	moodycamel::ConsumerToken ctok(candq); // Consumes from candq
	moodycamel::ProducerToken ptok(countq); // Produces for countq
	int mecount = 0;
	std::tuple<lint, int, lint, int>* retvals;
	bitset<tn> F;
	// Iterate until a kill-sentry is found
	while(true){
		candq.wait_dequeue(ctok, F); // Wait for new guy in queue

		if ((F.test(0) == 1) && (F.test(1) == 0)) { // Not a positive LTF
			std::ostringstream stream;
			stream << "Tester thread " << id << " terminating after testing " << mecount << " functions.\n";
			log(stream.str());
			return;		// Stop looking (termination signal)
		}

		// Holds the number of classes of each of 4 types associated with F
		retvals = new tuple<lint, int, lint, int>(0, 0, 0, 0);
		
		/* If an LTF, generate # of goldilocks functions in its orbit */
		/* Place into retvals */
		if (issep(F)) { 
			// Get the dual
			bitset<tn> Fd;
			dual(F, Fd); 

			// Store self-dualization as double-long array
			int FSD[2 * tn];
			for (int i = 0; i < tn; i++) {
				FSD[i] = F.test(i);
				FSD[i + tn] = Fd.test(i);
			}

			// And chow parameters
			// Since self-dual, the extra zeroth chow parameter is 2^(n+1)/2 = 2^n
			// Each other chow parameter of the self-dualization is double.
			int chow[n + 1];
			chowdualup(F, chow);		// The self-dualization has parameters double these
			
			// For each distinct anti-self-dualization
			bool newantisd = true;	// If this antiselfdual is distinct from those tested
			for (int i = 0; i <= n; i++) {
				if (!newantisd) {
					if (chow[i] != chow[i - 1]) {
						newantisd = true;
					}
				}

				if (newantisd) {
					// If self-dual pair
					if (chow[i] == tn / 2) {
						// Test xi = 0 for smallness
						bool isSmall = true;
						for (int j = 1; j <= n; j++) {	// testing the singleton values
							if ((FSD[two(j)]) && (j != i)) {
								isSmall = false;
								break;
							}
						}
						if (isSmall) {
							// Sn multiplicity computation
							// reps = size of Sn orbit of generator
							static lint max = fact(n);
							lint reps = max;
							int rchow[n];			// Get reduced chow parameters
							int p = 0;
							for (int j = 0; j <= n; j++) {
								if (j != i) {
									rchow[p] = chow[j];
									p++;
								}
							}

							int pcount = 1;
							for (int j = 1; j < n; j++) {
								if (rchow[j] == rchow[j - 1]) {
									pcount++;
								}
								else {
									reps /= fact(pcount);
									pcount = 1;
								}
							}
							reps /= fact(pcount);

							// Record number of classes from this generator
							get<0>(*retvals) += reps;
							get<1>(*retvals) += 1;
							get<2>(*retvals) += reps;
							get<3>(*retvals) += 1;
						}
					}
					else { // Not self dual (distinct)
						int numberPS = 2;			// Innocent until proven guilty

						// Test xi = 0 for smallness
						for (int j = 1; j <= n; j++) {	// testing the singleton values
							if ((FSD[two(j)]) && (j != i)) {
								numberPS--;
								break;
							}
						}

						// Test xi = 1 smallness
						int ti = two(i);
						for (int j = 1; j <= n; j++) {	// testing the singleton values
							if ((FSD[two(j) + ti]) && (j != i)) {
								numberPS--;
								break;
							}
						}

						// Now compute corresponding counts
						if (numberPS > 0) {
							
							// Sn multiplicity computation
							// reps = size of Sn orbit of the generator
							static lint max = fact(n);
							lint reps = max;
							int rchow[n];			// Get reduced chow parameters
							int p = 0;
							for (int j = 0; j <= n; j++) {
								if (j != i) {
									rchow[p] = chow[j];
									p++;
								}
							}

							int pcount = 1;
							for (int j = 1; j < n; j++) {
								if (rchow[j] == rchow[j - 1]) {
									pcount++;
								}
								else {
									reps /= fact(pcount);
									pcount = 1;
								}
							}
							reps /= fact(pcount);

							// Record number of classes from this generator
							if (numberPS == 2) {
								get<0>(*retvals) += reps;
								get<1>(*retvals) += 1;
								get<2>(*retvals) += reps;
								get<2>(*retvals) += reps;
								get<3>(*retvals) += 2;
							}
							else { // = 1
								get<2>(*retvals) += reps;
								get<3>(*retvals) += 1;
							}
						}
					}
					newantisd = false;
				}
			} // End for loop
		} // End F sep/ F testing and enumerating

		// Pass number of 
		countq.enqueue(ptok, *retvals);
		mecount++;

	} // End while loop
}

// Thread function: totals the counts resulting from the separate tests
void totaler(){
	lint tcount = 0;	        //Number of testcases.
	lint pcount = 0;			//Testcases*100
	lint GLcountSn = 0;						// Number of Hassett chambers quotiented by Sn
	lint GLcount = 0;
	lint PScountSn = 0;
	lint PScount = 0;

	int finishedcount = 0;
	
	moodycamel::ConsumerToken ctok(countq);
	log("Totaler: Initiated\n");
	
	while(tcount < TOTALT) {
		std::tuple<int, int, int, int> retvals;
		countq.wait_dequeue(ctok, retvals);
		
		// Tally the counts
		tcount++; pcount += 100;
		GLcount += std::get<0>(retvals);
		GLcountSn += std::get<1>(retvals);
		PScount += std::get<2>(retvals);
		PScountSn += std::get<3>(retvals);

		// Mark progress
		static lint percent = 1;
		if ((pcount/TOTALT) >= percent) {
			std::stringstream stream;
			stream << "Totaler: " << percent << "% complete.\n";
			stream << "Current progress:\n";
			stream << "n = " << n << "\n";
			stream << "Number Tested : " << tcount << "\n";
			stream << "Number Goldilocks(/Sn): " << GLcountSn << "\n";
			stream << "Number Goldilocks: " << GLcount << "\n";
			stream << "Number SemiGold(/Sn): " << PScountSn << "\n";
			stream << "Number SemiGold: " << PScount << "\n";
			output(stream.str());

			percent++;
		}
	}

	// Output results
	cout << "Final Results!" << endl;
	cout << "n = " << n << endl;
	cout << "Number Tested : " << tcount << endl;
	cout << "Number Goldilocks(/Sn): " << GLcountSn << endl;
	cout << "Number Goldilocks: " << GLcount << endl;
	cout << "Number SemiGold (/Sn): " << PScountSn << endl;
	cout << "Number SemiGold: " << PScount << endl;
	
	std::stringstream stream;
	stream << "Final Results!\n";
	stream << "n = " << n << "\n";
	stream << "Number Tested : " << tcount << "\n";
	stream << "Number Goldilocks(/Sn): " << GLcountSn << "\n";
	stream << "Number Goldilocks: " << GLcount << "\n";
	stream << "Number SemiGold(/Sn): " << PScountSn << "\n";
	stream << "Number SemiGold: " << PScount << "\n";
	output(stream.str());
}


// Main: original thread spawns others, and then reads functions into pool queue
int main() {
	// Real main begins here
	lessgreatinit(Great, Less);

	std::stringstream stream;
	stream << "Beginning execution at " << "\n";
	log(stream.str());

	// Initial thread produces for candq
	moodycamel::ProducerToken ptok(candq);

	// Creates an army of tester threads
	std::thread thdary[MAXTHREADS-2];
	for (int i = 0; i < MAXTHREADS-2; i++) {
		thdary[i] = std::thread(tester, i);

		std::stringstream stream;
		stream << "Main: spawned tester thread " << i << endl;
		log(stream.str());
	}

	// Creates the final thread which totals the various sums
	std::thread final = std::thread(totaler);
	log("Main: spawned totaler thread.\n");

	// Read each function from the file
	// Possibly buffered (TODO)
	ifstream infile;
	infile.open(readname, ios::in | ios::binary);

	if (!infile) {
		log("Read failure -- terminating.\n");
		return(1);
	}

	char* Fs = new char[tn+1];
	int i = 0;
	while (infile.read(Fs, tn+1)){
		Fs[tn] = '\0';
		bitset<tn> F = bitset<tn>(Fs, tn, '0', '1');

		while(candq.size_approx() > QUEUEMAX){ // Wait until queue has emptied
			std::this_thread::sleep_for(std::chrono::milliseconds(WAITFOR));
		} 
		candq.enqueue(ptok, F);
	}

	// Once all have been read, push MAXTHREADS-2 terminating tokens into candq
	for (int i = 0; i < MAXTHREADS-2; i++) {
		bitset<tn>* Fi = new bitset<tn>();
		Fi->set(0, true);

		candq.enqueue(ptok, *Fi); // Put the kill tokens in the queue
	}

	// Then wait on termination
	std::stringstream stream4;
	stream4 << "Main: " << std::to_string(MAXTHREADS-2);
	stream4 <<" kill-tokens sent, commencing wait.\n";
	log(stream4.str());
	for (int i = 0; i < MAXTHREADS-2; i++) {
		thdary[i].join();

		std::stringstream stream2;
		stream2 << "Main: received tester " << std::to_string(i);
		stream2 << " thread termination.\n";
		log(stream2.str());
	}

	final.join();
	log("Main: Terminating all execution.\n");
}