#include <Rcpp.h>
using namespace Rcpp;

// Author: Mathieu Lajoie, Date: September 2019

// Checks if both matrices have same dimensions
bool mat_dim_equal(SEXP mat1, SEXP mat2){
  NumericMatrix m1(mat1);
  NumericMatrix m2(mat2);
  if( (m1.ncol() != m2.ncol()) || (m1.nrow() != m2.nrow())){
    return false;
  }else{
    return true;
  }  
}

// [[Rcpp::export]]
IntegerMatrix cpp_obj_size(SEXP mseg){
  // input: seg matrix. 
  // output: seg matrix with object sizes
  
  // input
  IntegerMatrix xx(mseg);
  int nb_col = xx.ncol();
  int nb_row = xx.nrow();
  
  // output matrix
  IntegerMatrix res(nb_row, nb_col);
  // count vector
  IntegerVector counts(nb_row * nb_col);
  
  // loop through pixels and update counts
  for(int i=0; i< nb_row; i++){
    for(int j=0; j< nb_col; j++){
      if(xx(i,j) > 0){ // Positive pixels only?
        counts[xx(i,j)]++;
      }
    }
  }
  // fill output matrix
  for(int i=0; i< nb_row; i++){
    for(int j=0; j< nb_col; j++){
      res(i,j) = counts[xx(i,j)];
    }
  }
  
  return(res);
}

// [[Rcpp::export]]
NumericMatrix cpp_obj_mean(SEXP mseg, SEXP mval){
  // input: seg & value matrice 
  // output: seg matrix with object mean value
   
  // input
  IntegerMatrix ms(mseg);
  NumericMatrix mv(mval);
  
  int nb_col = ms.ncol();
  int nb_row = ms.nrow();
  
   if(!mat_dim_equal(mseg, mval)) {
     throw std::invalid_argument("Two matrices are not of the same dimension!");
   }
  
  // output matrix
  NumericMatrix res(nb_row, nb_col);
  
  // count vector
  IntegerVector counts(nb_row * nb_col);
  // value vector
  NumericVector values(nb_row * nb_col);
  
  // loop through pixels and update counts & values
  for(int i=0; i< nb_row; i++){
    for(int j=0; j< nb_col; j++){
      if(ms(i,j) > 0){ // Positive pixels only?
        counts[ms(i,j)] ++ ;
        values[ms(i,j)] += mv(i,j) ;
      }
    }
  }
  // fill output matrix
  for(int i=0; i< nb_row; i++){
    for(int j=0; j< nb_col; j++){
      if(ms(i,j) > 0){
        res(i,j) = values[ms(i,j)]/counts[ms(i,j)];
      }
    }
  }
  
  return(res);
}

// [[Rcpp::export]]
IntegerMatrix cpp_obj_contour(SEXP mseg, int neisize = 8){
  // input: seg matrix
  // output: seg matrix restricted to (inside) contour 

  if(neisize!=4 && neisize!=8){
    throw std::invalid_argument("neisize must be either 4 or 8");
  }
  
  //Nei type:  .....4.....    .....8.....
  int NI[8] = { 0, 0,-1, 1,  -1,-1, 1, 1};
  int NJ[8] = {-1, 1, 0, 0,  -1, 1,-1, 1};
  
  // input matrix
  IntegerMatrix xx(mseg);
  int nb_col = xx.ncol();
  int nb_row = xx.nrow();
  
  // output matrix
  IntegerMatrix res(nb_row, nb_col);
  
  // for each pixel
  for(int i=0; i< nb_row; i++){
    for(int j=0; j< nb_col; j++){
      
      int pix_val = xx(i,j);
      
      if(pix_val > 0){ // Positive pixels only
        
        bool pix_contour = false;
        
        // for each nei (4 or 8)
        for(int x = 0; x < neisize; x++){
          int ni = i + NI[x];
          int nj = j + NJ[x];
          if( ni >= 0 && ni < nb_row && nj >= 0 && nj < nb_col){
            if(xx(ni,nj) !=  pix_val){
              pix_contour = true;
            }
          }
        }
        if(pix_contour){
          res(i,j) = pix_val;
        }
      }
    }
  }
  return(res);
}

