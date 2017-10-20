// GoldilocksEnumParallel.cpp
// by Nicolle Gruzling 
// Adapted by Connor Halleck-Dube as part of SUMRY 2017
// https://arxiv.org/abs/1709.03663

// This implements the first part of an algorithm for the enumeration of 
// Goldilocks linear threshold functions (GLTFs). The program generates all 
// hypercomplete boolean functions on n variables; writes them to file.
// These are then read and tested for separability by GoldilocksTestParallel.cpp

// This piece of the algorithm can be found in:
// R. O. Winder. Enumeration of seven-argument threshold functions. 
// 		IEEE Transactions on Electronic Computers, EC-14(3):315–325, 1965.


#include "usefcns.h"
#include "bigint.h"
#include "blockingconcurrentqueue.h"
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


//const unsigned n=3; const unsigned tn = 8;
//const unsigned n=4; const unsigned tn = 16;
//const unsigned n=5; const unsigned tn = 32; 
//const unsigned n=6; const unsigned tn = 64; 
//const unsigned n=7; const unsigned tn = 128; 
//const unsigned n=8; const unsigned tn = 256; 
//
const unsigned n=9; const unsigned tn = 512;

vector<int> Great[tn],Less[tn];

#include "functions.cpp"

// Store output 
char outname[] = "/home/fas/payne_sam/cjh69/project/GoldCands9.dat";

int buflen; 						// Write buffer length (in bitsets)
const int bufsize = 2097152; 		// Write buffer size (in chars) 
									// 		Must be mult of tn, 
									//		smaller version 65536
int bufi = 0; 						// Number of bitsets currently in buffer

char buffer[bufsize];				// Write buffer

// Buffered write of bitset data to file
void write(ofstream& outfile, bitset<tn>& F) {
	static int buflen = bufsize*8/tn;
	static const int charperbs = tn/8;
	for (int ci = 1; ci <= charperbs; ci++){
		char c = 0;
		for (int bi = 7; bi >= 0; bi--) {
			if (F.test(tn - ci*8 + bi)){
				c |= 0x1;
			}
			c <<= 1;
		}

		buffer[bufi*charperbs + ci - 1] = c;
	}
	bufi++;

	if (bufi == buflen) { // Flush buffer
		outfile.write(buffer, bufsize);
		bufi = 0;
	}
}

// Final flush of buffer
void flush (ofstream& outfile) {
	outfile.write(buffer, bufi);
}

void initless(bitset<tn> less[]);   //Initializes in.

int main(){
	ofstream outfile(outname, ios::binary);

	lessgreatinit(Great,Less);
	bitset<tn> lessa[tn]; initless(lessa);

	lint tcount=0,				        	//Number of testcases.
	BigInt septotal(tcount);

	vector< bitset<tn> > stk;           	//A stack to hold fcns to process
	bitset<tn> F,free; 
	F.reset(); free.set();

	stk.push_back(F);  						//Push the first fcn on the stack. 
	stk.push_back(free);					//Along with its free elements.

 	while(!stk.empty() ){
		free = stk.back();stk.pop_back();	//Pop the top set and it's free 
		F = stk.back(); stk.pop_back(); 	//variables off the stack.
		while(free.count()>0){      
			int j = tn-1;         			//Find the largest free element.
			while(!free.test(j))
				j--;
			stk.push_back(F);         		//Push a copy of F onto the stack.
			stk.push_back(free);

			free &= lessa[j];         		//Remove elements less than j.

		}
		tcount++;
		write(outfile, F);

		if(!stk.empty()){
			free = stk.back();stk.pop_back(); //Pop the top set and it's free 
			F = stk.back(); stk.pop_back();   //variables off the stack.

			unsigned j = tn-1;        		//Find the largest free element.
			while(!free.test(j))
				j--;
			F.set(j); free.reset(j);
			if(posn(j,n-1) && posn(j,n-2) ){
				unsigned z = comp(n-2,j); set(z,n-1);
				free&=lessa[z];
			}
			stk.push_back( F );				//Push the first set on the stack. 
			stk.push_back( free );        	//Along with it's free elements.
		}
	}

	flush(outfile);
	outfile.close();

	cout<<"\nNumber Generated : "<<tcount<<endl;
	return 0;
}

//*************************************************************************
//*************************************************************************

void initless(bitset<tn> less[]){      
	for(int i=0;i<tn;i++){
		less[i]=bitset<tn>();
		less[i].set();
		less[i].reset(i);
		for(int j=0;j<tn;j++){
			if( lessdot(n,j,i) )
				less[i].reset(j);
		}
	}
}

