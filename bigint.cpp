#include "bigint.h"
BigInt::BigInt( unsigned long long val){
  size=30;
  unsigned max=~0; 
  base=1;
  while(base*10<max)
    base*=10;

  for(int i=0;i<size;i++)
    number[i]=0;
    
  for(int i=0;val!=0 && i<size; i++){
    number[i] = val%base;
    val/=base;
  } 
}

BigInt::BigInt( const BigInt& rhs ){
  size=rhs.size;
  base=rhs.base;
  for(int i=0;i<size;i++)
    number[i]=rhs.number[i];
}

BigInt::BigInt(int arr[]){
  size=30;
  unsigned max=~0; 
  base=1;
  while(base*10<max)
    base*=10;
  for(int i=0;i<size;i++)
    number[i]=arr[i];
  
}
void BigInt::getArray( int arr[]){
  for(int i=0;i<size;i++)
    arr[i]=number[i];
}

const BigInt &BigInt::operator=(const BigInt & rhs){
    for( int i= 0; i<size; i++ )
      number[i] = rhs.number[i];
    return *this;  
}
  
BigInt BigInt::operator+=( BigInt& rhs){
  unsigned long long carry=0;
  for(int i=0;i<size;i++){
    unsigned long long sum = number[i]+rhs.number[i]+carry;
    if(sum >= base){
      unsigned long long rem = sum%base; 
      number[i] = rem;
      carry=sum/base;
    }
    else{
      number[i]+= rhs.number[i]+carry;
      carry=0;
    }  
  }   
  return *this;
}
BigInt BigInt::operator+=( unsigned long long rhs){
  BigInt t(rhs);
  (*this)+=t;
  return *this;
} 

bool BigInt::operator==( const BigInt & r) const{
  for(int i=0;i<size;i++)
    if(number[i]!=r.number[i]) 
      return false;
  return true;
}

ostream &operator<<(ostream& out, BigInt& rhs){
  int i=rhs.size-1;
  while(i>0 && rhs.number[i]==0)
    i--;
  out<<rhs.number[i--];
  
  unsigned max=rhs.base-1;
  int num_digits=0;
  while(max/10 > 0 ){
    num_digits++;
    max/=10;
  }  
  for(;i>=0;i--){
    out.setf( ios::right, ios::adjustfield );
    out<<setw(num_digits+1)<<setfill('0')<<rhs.number[i];
  }  
  return out;
}
