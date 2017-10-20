// functions.cpp
// Written by Nicolle Gruzling
// Comments by Nicolle Gruzling and Connor Halleck-Dube
// Contains an assortment of functions for computation on threshold functions
// More information: 
// http://unbc.arcabc.ca/islandora/object/unbc%3A15882/datastream/PDF/download/citation.pdf

// Boolean functions F: {0, 1}^n --> {0, 1} encoded as bitstrings of length 2^n
// Given F, the entry F[i] encodes the value of F on the binary repn of i.
// E.g. if n=6, then F[4] = 1 means that F(000100) = 1.

using namespace std;

// Compute various Chow parameters of boolean function F
void chowa(const bitset<tn>& F, int a[]){
  for(int i=0;i<n;i++)
    a[i]=0;
    
  for(int i=0;i<tn;i++){
    if(F.test(i)){
      for(int j=0;j<n;j++){
        a[j]+=posn(i,j);
      }
    }
  }  
  return;
}
void chowav(const bitset<tn>& F, vector<int>& a){
  for(int i=0;i<n;i++)
    a.push_back(0);
    
  for(int i=0;i<tn;i++){
    if(F.test(i)){
      for(int j=0;j<n;j++){
        a[j]+=posn(i,j);
      }
    }
  }  
  return;
}
void chowdualup(const bitset<tn>& F, int a[]){
  for(int i=0;i<=n;i++)
    a[i]=0;
    
  for(int i=0;i<tn;i++){
    if(F.test(i)){
      for(int j=0;j<n;j++){
        a[j]+=posn(i,j);
      }
    }
    else{
      int h=comp(n+1,i);
      for(int j=0;j<=n;j++){
        a[j]+=posn(h,j);
      }
    }
  }  
  return;
}

// Compute the number of boolean functions which can be reached from F
// by permutation and complementation of any arguments of F
lint reps(const bitset<tn>& F, int m){
  static lint max = fact(n)*tn;
  lint rep = max;
  int a[n]; chowa(F,a);
  bool even = m%2==0;
  int pcount = 1;
  rep = ((even && a[n-1]==m/2)? rep/2 : rep);
  
  for(int i=n-2;i>=0;i--){
    rep = ((even && a[i]==m/2)? rep/2 : rep);
    
    if(a[i]==a[i+1]){
      pcount++;
    }  
    else{
      rep/=fact(pcount);
      pcount=1; 
    }   
  }
  rep/=fact(pcount);
  return rep;      
}
// Compute the number of boolean functions which can be reached from F
// by permutation, complementation, or self-dualization or anti-self-dualization
lint repsdualup(const bitset<tn>& F){
  int m=tn;
  static lint max = fact(n+1)*tn*2;
  lint rep = max;
  int a[n+1]; chowdualup(F,a);
  int pcount = 1;
  rep = ((a[n]==m/2)? rep/2 : rep);
  
  for(int i=n-1;i>=0;i--){
    rep = ((a[i]==m/2)? rep/2 : rep);
    
    if(a[i]==a[i+1]){
      pcount++;
    }  
    else{
      rep/=fact(pcount);
      pcount=1; 
    }   
  }
  rep/=fact(pcount);
  return rep;      
}

//Returns true if i (assumed in F) is a border point of F.
bool isinborder(unsigned i,const bitset<tn>& F){ 
  for(int j=0;j<n;j++){                     
    unsigned less=i;
    if(posn(i,j)==1){
      set(less,j,0);
      if( F.test(less) )
        return false;
    }  
  }    
  return true;
}
//Returns true if i (assumed not in F) is a border point of F.
bool isoutborder(unsigned i,const bitset<tn>& F){ 
  for(int j=0;j<n;j++){                        
    unsigned great=i;
    if(posn(i,j)==0){
      set(great,j,1);
      if( !F.test(great) )
        return false;
    }
  }    
  return true;
}

//Returns true if i (assumed in F) is a boundary point of F.
bool ishighbound(int i, const bitset<tn>& F){ 
  for(int j=0;j<less[i].size();j++)           
    if( F.test( less[i][j] ) )
      return false;
  return true;
}
//Returns true if i (assumed not in F) is a boundary point of F.
bool islowbound(int i, const bitset<tn>& F){  
  for(int j=0;j<great[i].size();j++)         
    if( !F.test( great[i][j] ) )
      return false;
  return true;
}

// Tests whether F is 2-monotonic
bool ismonotonic(int an,int atn,const bitset<tn>& b){
  bool val = 1;                               //In most cases we always testing
                                              //sets that that have less 1's
					      //than 0's.  
  unsigned mask = 0;                          //We never care about the bits 
  for(int i=0;i<an;i++)                       //greater than n, mask is used
    set(mask,i,1);                            //make sure they are 0.
    
  unsigned h=an/2;  				
    
  for(unsigned m=0; m<atn; m++){              //m is used to find the different
                                              //faces of the hypercube.
					   
    if(count(an,m,1)<=2 && count(an,m,1)<=h){
    				      
      for(unsigned i=0; i<atn; i++){
        int icomp = (~(i^m))&mask;
        if(b.test(i)==val && b.test(icomp)==val){//In the m face, if i and it's
                                              //pair are in the set, then check
					      //if another pair is not in the
					      //set.
					      
	  for(unsigned k=0; k<icomp; k++ ){
	    if( ((k&m)==(i&m)) &&             //Check k is in same face as i.
	        b.test(k)!=val &&             //Check that k is in opposite set.
	        b.test((~(k^m))&mask)!=val ){ //Check that the pair of k is in
	                                      //the same set as k.
	      return false;                   //If so b is not monotonic.
	    }
	  }
        } 
      }
    }
  }
  return true;                                //Innocent until proven guilty.
}
bool lexlesseq(double* b1,double* b2,int s){
  int i=0;
  while(i<s && (b1[i] - b2[i] == 0) )
    i++;
  if(i==s) 
    return true;
  if( b1[i]-b2[i] < 0)
    return true;
  return false;     
}

// Tests a linear inequality system mat for solution by the simplex method
// True if a solution exists, false otherwise
bool dual_simplex(double** mat,const int p,const int q, double* soln){
do{ 
  static double e=0.000000001; 
  bool opt = true;
  int i=0;
  for(;i<p+q;i++){
    if(mat[i][0]<0){
      opt = false;
      break;
    }
  }
  if(opt){
    soln[0] = mat[0][0];
    for(int k=1;k<q;k++){
      soln[k] = mat[p+k][0];
    }
    return true;
  }
  else{
    bool issoln = false;
    for(int j=1;j<q;j++){
      if(mat[i][j]>e){
        issoln = true;
        break;
      }
    }
    if(!issoln)
      return false;
    
    vector<double*> bs; vector<int> bp;
    for(int j=1; j<q; j++){
      if(mat[i][j]>e){
        double* b = new double[p+q];
	for(int k=0;k<p+q;k++){
	  b[k]=mat[k][j]/mat[i][j];
	}  
	bs.push_back(b);
	bp.push_back(j);
      }
    }
    int j;
    for(int k=0;k<bs.size();k++){
      bool isless = true;
      for(int l=0;l<bs.size();l++){
        if(!lexlesseq(bs[k],bs[l],p+q)){
	  isless = false;
	}  
      }
      if(isless){
        j=bp[k];
	for(int m=0;m<p+q;m++){
	  mat[m][j]=bs[k][m];
	}
        break;
      }
    }
    
    for(int k=0;k<bs.size();k++)
      delete bs[k];
    
    double b[p+q];
    for(int k=0;k<p+q;k++)
      b[k]=mat[i][k];
    
    for(int l=0;l<q;l++)
      if(l!=j)
        for(int k=0;k<p+q;k++)
          mat[k][l]= mat[k][l] - b[l]*mat[k][j];
    
  }
}while(true);  
}

// Initializes the "Less" and "Great" arrays
// They encode a partial order relationship between vectors in {0, 1}^n
// R. O. Winder. Enumeration of seven-argument threshold functions. 
//       IEEE Transactions on Electronic Computers, EC-14(3):315–325, 1965.
void lessgreatinit(vector<int> great[], vector<int> less[]){	
  for(int i=0;i<tn;i++){
    for(int j=i+1;j<tn;j++){
      if(  lessdot(n,i,j) ){
        bool tightest = true;
        for(int k=0;k<great[i].size();k++){
	  if( lessdot( n,(great[i])[k],j) ){
	    tightest=false;
	  }
	}
	if(tightest){
	  great[i].push_back(j);
	}
      }
    }
    for(int j=i-1;j>=0;j--){
      if( lessdot(n,j,i) ){
        bool tightest = true;
        for(int k=0;k<less[i].size();k++){
	  if( lessdot( n,j,less[i][k]) ){
	    tightest=false;
	  }
	}
	if(tightest){
	  less[i].push_back(j);
	}
      }
    }
  }    
}

// Tests whether a boolean function F is a linear threshold function
bool issep(bitset<tn>& F){
  vector< vector<int> > constraints;
  for(int i=0;i<tn;i++){
    if(F.test(i)){
      if(ishighbound(i,F)){
        vector<int> a;
	a.push_back(0);
	a.push_back(-1);
        for(int j=0; j<n; j++)
          a.push_back( static_cast<int>(posn(i,j)) );
        constraints.push_back(a);
      }
    }
    else{
      if(islowbound(i,F)){
        vector<int> a;
	a.push_back(-1);
	a.push_back(1);
        for(int j=0; j<n; j++)
          a.push_back( -1*static_cast<int>(posn(i,j)) );
        constraints.push_back(a);
      }
    }
  }
  
  int num_cols = n+2;  
  int num_rows = constraints.size()+ n-1; 	
  
  double **mat = new double*[num_rows+num_cols];
  
  int i=1;
  mat[0] = new double[num_cols];
  mat[0][0]=0;
  for(int j=1;j<num_cols;j++)
    mat[0][j]=1.0;
    
  for(;i<=constraints.size();i++){
    mat[i] = new double[num_cols];
    for(int j=0;j<num_cols;j++){
      mat[i][j]= static_cast<double>(constraints[i-1][j]);
    }
  }
  for(int p=2;i<=num_rows;i++,p++){
    mat[i] = new double[num_cols];
    for(int j=0;j<num_cols;j++){
      mat[i][j]=(j==p?-1:(j==p+1?1:0));
    }  
  }
  for(int p=1;i<num_rows+num_cols;i++,p++){
    mat[i] = new double[num_cols];
    for(int j=0;j<num_cols;j++)
      mat[i][j]=(j==p?1:0);
  }  
  double soln[num_cols];
  bool sep = dual_simplex(mat,num_rows,num_cols,soln);

  for(i=0;i<num_rows+num_cols;i++)
    delete mat[i];
  delete mat; 
  
  return sep;
}
bool issep(bitset<tn>& F, double soln[]){
  vector< vector<int> > constraints;
  for(int i=0;i<tn;i++){
    if(F.test(i)){
      if(ishighbound(i,F)){
        vector<int> a;
	a.push_back(0);
	a.push_back(-1);
        for(int j=0; j<n; j++)
          a.push_back( static_cast<int>(posn(i,j)) );
        constraints.push_back(a);
      }
    }
    else{
      if(islowbound(i,F)){
        vector<int> a;
	a.push_back(-1);
	a.push_back(1);
        for(int j=0; j<n; j++)
          a.push_back( -1*static_cast<int>(posn(i,j)) );
        constraints.push_back(a);
      }
    }
  }
  
  int num_cols = n+2;  
  int num_rows = constraints.size()+ n-1; 	
  
  double **mat = new double*[num_rows+num_cols];
  
  int i=1;
  mat[0] = new double[num_cols];
  mat[0][0]=0;
  for(int j=1;j<num_cols;j++)
    mat[0][j]=1.0;
    
  for(;i<=constraints.size();i++){
    mat[i] = new double[num_cols];
    for(int j=0;j<num_cols;j++){
      mat[i][j]= static_cast<double>(constraints[i-1][j]);
    }
  }
  for(int p=2;i<=num_rows;i++,p++){
    mat[i] = new double[num_cols];
    for(int j=0;j<num_cols;j++){
      mat[i][j]=(j==p?-1:(j==p+1?1:0));
    }  
  }
  for(int p=1;i<num_rows+num_cols;i++,p++){
    mat[i] = new double[num_cols];
    for(int j=0;j<num_cols;j++)
      mat[i][j]=(j==p?1:0);
  }  
  bool sep = dual_simplex(mat,num_rows,num_cols,soln);

  for(i=0;i<num_rows+num_cols;i++)
    delete mat[i];
  delete mat; 
  
  return sep;
}
