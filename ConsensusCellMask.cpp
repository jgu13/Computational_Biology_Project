#include <Rcpp.h>
//#include <math.h>
  using namespace Rcpp;

struct obj_pixel{
  int i;
  int j;
  int obj_id; //object id
  double prob; //object prob
  char marker;
};

bool compare_obj_pixel(const obj_pixel &a, const obj_pixel &b){
  if(a.prob == b.prob){
    return a.obj_id > b.obj_id; //give preference to obj with higher id
  }else{
    return a.prob > b.prob;
  }
}

// fill pixel with object probability mean
void average_obj_prob(NumericMatrix s, NumericMatrix p){
  int nrow = s.nrow(), ncol=s.ncol();
  
  double obj_prob_sum[nrow*ncol];
  int pix_count[nrow*ncol];
  
  // loop through every pixel of the probability map
  for(int i=0; i < nrow; i++){
    for(int j=0; j < ncol; j++){
      int id = s(i,j);
      if(id > 0){
        obj_prob_sum[id]+= p(i,j);
        pix_count[id]++;
      }
    }
  }
  
  for(int i=0; i<nrow; i++){
    for(int j=0; j<ncol; j++){
      int id = s(i,j);
      double obj_mean = obj_prob_sum[id]/pix_count[id];
      if(id > 0){
        p(i,j) = obj_mean;
      }
    }
  }
}

// [[Rcpp::export]]
NumericMatrix cpp_segmentation_merge(SEXP segA, SEXP probA, SEXP segB, SEXP probB){
  
  NumericMatrix mSA(segA);
  NumericMatrix mPA(probA);
  NumericMatrix mSB(segB);
  NumericMatrix mPB(probB);
  
  // TODO:
    // Check that all matrices have the same dimensions (going outside bounds will crash!)
  if( (mSA.ncol() != mPA.ncol()) || (mSA.nrow() != mPA.nrow())){
    throw std::invalid_argument("Two matrices are not of the same dimension!");
  } 
  
  int nb_row = mSA.nrow(), nb_col=mSA.ncol();
  
  // Init result matrix
  NumericMatrix mRES(nb_row, nb_col);
  
  // Init vector of pixels
  std::vector<obj_pixel> vPix;
  
  //average probability
  average_obj_prob(mSA, mPA);
  average_obj_prob(mSB, mPB);
  
  // Add pixels to vector   
  for(int i=0; i< nb_row; i++){
    for(int j=0; j< nb_col; j++){
      
      obj_pixel pa = obj_pixel();
      pa.prob = mPA(i,j);
      pa.obj_id = mSA(i,j);
      pa.i = i; pa.j = j;
      pa.marker='A';
      vPix.push_back(pa);
      
      obj_pixel pb = obj_pixel();
      pb.prob = mPB(i,j);
      pb.obj_id = mSB(i,j);
      pb.i = i; pb.j = j;
      pb.marker='B';
      vPix.push_back(pb);
    }
  }
  
  //Rcout<<"init vector done"<<"\n";
  std::sort(vPix.begin(), vPix.end(), compare_obj_pixel);
  
  //loop through pixels
  int k = 0;
  
  while(k < vPix.size()){
    obj_pixel& p = vPix[k];
    
    int obj_id = p.obj_id; //current obj_id
    bool conflict = false;
    
    // check conflict for current object
    int k2 = k; 
    while(k2 < vPix.size() && vPix[k2].obj_id == obj_id){
      
      if(mRES(vPix[k2].i , vPix[k2].j) != 0){
        conflict = true;
      }
      
      k2++;
    }
    
    // if no conflict, we add this object
    if(!conflict){
      while(k < vPix.size() && k < k2){
        mRES(vPix[k].i, vPix[k].j) = (vPix[k].marker == 'A')? 2:3;
        k++;
      }
    }
    else{
      k = k2; //skip this object 
    }
  }
  
  return(mRES);
}












