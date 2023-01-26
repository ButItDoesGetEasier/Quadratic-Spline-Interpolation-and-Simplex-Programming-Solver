#source("DamianEx03.R") # import to get access of AugCoeffMatrix()

pivot <- function(augMatCoeff, i, j) { # swaps the ith and jth row of augmatcoeff
  temp = augMatCoeff[i,] # store ith row in temp
  augMatCoeff[i,] = augMatCoeff[j,] # put jth row in the ith
  augMatCoeff[j,] = temp # put ith row in the jth
  return(augMatCoeff)
}

GaussianElimination <- function(augMat) {
  augMatCoeff = augMat$augcoeffmatrix
  dne = FALSE
  for (i in 1:(length(augMat$variables)-1)) { # if n=num of vars, there are n-1 iterations
    largest = i # find biggest element in pivot column, set initially at i
    for (j in (i+1):(length(augMat$variables))) { # loop through every row elements below i
      if (abs(augMatCoeff[j,i]) > abs(augMatCoeff[largest, i])) { # check if absolute val of checked element is larger than saved element in largest
        largest = j # if so, re-set largest to j
      }
    }
    
    if (augMatCoeff[largest, i] == 0) { # if the largest is 0, there are no solution
      dne = TRUE;
      break
    } 
    augMatCoeff = pivot(augMatCoeff, i, largest) # swap ith row to the row containing the largest jth elem
    
    for (k in (i+1):(length(augMat$variables))) { # loop from i+1 to however many vars there are
      pivotElem = augMatCoeff[i,i] # the [i,i]th elem is the pivot element
      multiplier = augMatCoeff[k,i]/pivotElem # multiplier is the element [k,i] divided by the pivot elem
      temp = augMatCoeff[i,] * multiplier # set temp vector whose val is the pivot row multiplied by the multiplier
      augMatCoeff[k,] = augMatCoeff[k,]-temp # subtract temp from the kth row
    }
  }
  if(dne) { # if there are no sol'n, set NA as the solution set
    result = list(solutionSet=NA,Variables=augMat$variables)
  } else {
    #backward substitution
    x = numeric(length(augMat$variables)) # vector to contain values of unknown vars
    for (i in (length(augMat$variables)):1) { # loop from however many vars there are to 1 (top to bottom)
      if (i==length(augMat$variables)) { # if i is the number of vars, 
       # ith =                rhs                         - sum of coefficient*solved variable value pairs                                   / diagonal element of i
        x[i] = (augMatCoeff[i,length(augMat$variables)+1] - sum(augMatCoeff[i, i:length(augMat$variables)] * x[i:length(augMat$variables)])) / augMatCoeff[i,i]
      } else {
        x[i] = (augMatCoeff[i,length(augMat$variables)+1] - sum(augMatCoeff[i, (i+1):length(augMat$variables)] * x[(i+1):length(augMat$variables)])) / augMatCoeff[i,i]
      } # i separated them in if-cases since i starts at the biggest index, which would make x[(i+1):length(augMat$variables)] go out of index
        # either way, there is no need to reach that index if i is the biggest
    }
    result = list(solutionSet=x,Variables=augMat$variables, matrix=augMatCoeff) # return labelled list containing the solution set, vars, and the used matrix
   }
  return(result)
}

GaussJordanElimination <- function(augMat){
  augMatCoeff = augMat$augcoeffmatrix
  dne = FALSE
  for(i in 1:length(augMat$variables)) { # loop from 1 to however many vars/eqs there are (since every eq would be a pivot row)
    if (i != length(augMat$variables)) { # pivoting would not happen for the last iteration
      largest = i # finds biggest element in column, pivots; works the same way as the one from gaussian elimination
      for (j in (i+1):(length(augMat$variables))) { # find diag dominant element on column i
        if (abs(augMatCoeff[j,i]) > abs(augMatCoeff[largest, i])) {
          largest = j
        }
        
      }
      if (augMatCoeff[largest, i] == 0) { # works the same as gaussian elimination
        dne = TRUE;
      } 
      augMatCoeff = pivot(augMatCoeff, i, largest)
    }
    if(dne) {
      result = list(solutionSet=NA,Variables=augMat$variables)
    } else {
      augMatCoeff[i,] = augMatCoeff[i,] / augMatCoeff[i,i] # normalize ith 
      for (j in 1:length(augMat$variables)) { # loop from 1 to however many vars there are
        if (j == i) {next} # if i = j, skip current iteration (pivot row would not modify itself)
        else {
          multiplier = augMatCoeff[j,i] # multiplier is [j,i]th element 
          temp = augMatCoeff[i,] * multiplier # set temp to contain the ith row multiplied by the multiplier
          augMatCoeff[j,] = augMatCoeff[j,] - temp # modify jth row: subtract temp from it
        }
      }
      # return labelled list containing sol'n set (the rhs of each rows), vars, and the used matrix
      result = list(solutionSet=as.vector(augMatCoeff[1:length(augMat$variables),length(augMat$variables)+1]),Variables=augMat$variables, matrix=augMatCoeff)
    }
  }
  return(result)
}


# E1 <- function (x1, x2, x3) 0.1 * x1 + 7 * x2 + -0.3 * x3 + 19.3 # test case
# E2 <- function (x1, x2, x3) 3 * x1 + -0.1 * x2 + -0.2 * x3 + -7.85
# E3 <- function (x1, x2, x3) 0.3 * x1 + -0.2 * x2 + 10 * x3 + -71.4
# system <- list(E1, E2, E3)
# augMat = AugCoeffMatrix(system)
# print(augMat)
# 
# print("GAUSSIAN ELIMINATION")
# print(GaussianElimination(augMat))
# print("GAUSS JORDAN ELIMINATION")
# print(GaussJordanElimination(augMat))