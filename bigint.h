#ifndef BIGINT_H
#define BIGINT_H

#include "stdafx.h"
#include <iostream>
#include <iomanip>

using namespace std;

class BigInt{
  friend ostream &operator<<(ostream&,BigInt&);
public: 
  BigInt( unsigned long long = 0);
  BigInt( const BigInt& rhs );
  BigInt( int[] );
  void getArray( int[] );
  const BigInt &operator=(const BigInt &);
  BigInt operator+=( BigInt&);
  BigInt operator+=( unsigned long long );
  bool operator==( const BigInt &) const;
  bool operator!=( const BigInt & r ) const{ return !(*this==r); }   
private:  
  unsigned number[30];
  int size;
  unsigned long long base;
};

#endif
