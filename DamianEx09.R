# assume that values in the tableau are of one type of constraint (i.e. only <= or only >=). Maximization: only <=. Minimization: only >=.
# if there are mixed constraints in the original problem, values in unmatched constraints are multiplied by -1 to flip inequality.
simplex <- function(tableau, isMax, problem) { 
  # all constraint rows are assumed to be <=, if they're >=, values were multiplied by -1 first before being put in the tableau
  if(isMax==FALSE) { # if minimization, tableau has to be modified
    newMat = matrix( # create matrix with same num of rows as tableau but with no columns yet 
      nrow=nrow(tableau),
      ncol=0,
    ) # the plan is to remove slack and z columns, transpose the x values and Ans columns, then add the slacks and z back
    labels = colnames(tableau) # get colnames of tableau
    xAndZ = 0 # will signify number of x values + z in tableau
    for(i in 1:length(labels)) { # loop through colnames
      if (substring(labels[i],1,1)=="x" || i == length(labels)) { # if it starts with x or is the last colname (i.e. Ans)
        newMat = cbind(newMat, tableau[,i]) # add the ith column to the newMat
        xAndZ = xAndZ + 1 # increment xAndZ
      }
    }
    newMat[,ncol(newMat)] = newMat[,ncol(newMat)] * -1 # multiply last column and row by -1
    newMat[nrow(newMat),] = newMat[nrow(newMat),] * -1
    newMat[nrow(newMat),ncol(newMat)]=0 # set value in the farthest low-east to 0
    newMat = t(newMat) # transpose the matrix
    
    slacksMat <- matrix( # matrix to store the slack columns
      data=0, # initially all 0
      nrow=nrow(newMat), # same nrow as newMat
      ncol=xAndZ # since it's minimization, the number of slack columns (plus z) relies on xAndZ
    )
    for(i in 1:ncol(slacksMat)) { # set i,ith elements in slacksMat to 1
      slacksMat[i,i] = 1
    }
    
    colLabel = c() # set up the column labels here (S1, S2, ... Sn, x1, .. xn, Z, Ans)
    for(i in 1:(ncol(tableau)-xAndZ-1)) { # perform as many times as number of cols minus xAndZ minus 1 (for Ans)
      colLabel = append(colLabel, paste("S",i,sep="")) # append Si
    }
    for(i in 1:xAndZ) { # append x and z values
      if (i == xAndZ) { # last iteration is for z
        colLabel = append(colLabel, paste("Z",sep=""))
      } else { # else, xi
        colLabel = append(colLabel, paste("x",i,sep=""))
      }
    }
    colLabel = append(colLabel, paste("Ans",sep="")) # add Ans
    # the newMat is composed of x values, followed by slacks, then the ans column
    newMat = cbind(newMat[,1:(ncol(newMat)-1)],slacksMat,newMat[,ncol(newMat)]) 
    tableau=newMat # set newMat as tableau
    colnames(tableau) = colLabel # set colLabel as tableau's colnames
    rownames(tableau) = 1:nrow(tableau) # set rownames as a counting from 1 up
  }
  while(TRUE){ # perform while there are negative elements in the last row
    pivotCol = -1 # set both to -1 to signify being invalid
    pivotRow = -1
      
    for (i in 1:(ncol(tableau)-1)) { # loop through every column except the last (the last col is permitted to have negative nums)
      # if the last value of the ith column is negative AND (pivotCol is unmodified OR the ith value is less than the pivotColth)
      if ((tableau[nrow(tableau), i] < 0) && ((pivotCol == -1)||(tableau[nrow(tableau), i] < tableau[nrow(tableau), pivotCol]))){
        pivotCol = i # set i as pivotCol
      }
    }
    if(pivotCol == -1) {break} # if pivotCol is unmodified, there're no more negative values in the last row. Simplex is finished.
      
    for (i in 1:(nrow(tableau)-1)) { # loop through every value in the pivotCol
      # if the value in the last column and pivotCol are non-negative AND (pivotRow is unmodified OR the ith test ratio is less than the stored one)
      if ((tableau[i,ncol(tableau)] > 0 && tableau[i,pivotCol] > 0) && ((pivotRow == -1)||
          (tableau[i,ncol(tableau)]/tableau[i,pivotCol] < tableau[pivotRow,ncol(tableau)]/tableau[pivotRow,pivotCol]))) {
        pivotRow = i # set i as pivotRow
      }
    }
    if(pivotRow == -1) { # if pivotRow is unmodified, there are no valid test ratios. Hence, no feasible solution
      # if problem is true, add a shipping.num label. All but the final tableau states No Feasible Solution
      if(problem) {list <- list(final.tableau = tableau, basic.solution = "No Feasible Solution", opt.val="No Feasible Solution",shipping.num="No Feasible Solution")}
      else {list <- list(final.tableau = tableau, basic.solution = "No Feasible Solution", opt.val="No Feasible Solution")}
      return(list) # return list
    }
    tableau[pivotRow,] = tableau[pivotRow,]/tableau[pivotRow,pivotCol] # normalize the pivotRow
    for (i in 1:nrow(tableau)) { # loop through every row
      if (i == pivotRow) {next} # skip pivotRow
      else {
        # set ith row as the ith row - (pivotRow * the element in ith row aligned with the pivotCol)
        tableau[i,] = tableau[i,] - (tableau[pivotRow,] * tableau[i,pivotCol]) 
      }
    }
  }
  
  soln <- c() # store final solution here
  if(isMax) { # if it's a maximization problem, solution depends on the answer column
    for (i in 1:(ncol(tableau)-1)) { # loop through every column except the last
      oneth = -1 # initially -1, will signify the row number of the single 1 in a column of zeroes. will be -1 if not cleared
      for (j in 1:nrow(tableau)) { # loop through every value in the column (nested loop)
        # if (the j,ith value is not 0 or 1) or (the j,ith value is not 0, can be 1, but the oneth is already set. (possibly two 1's in a column))
        if ((tableau[j,i] != 0) && ((tableau[j,i] != 1)||(oneth != -1))) { 
          oneth = -1 # set oneth as -1 then break
          break
        } else if (tableau[j,i] == 1) { # if it's a lone (so far) 1, set j as oneth
          oneth = j
        }
      }
      if (oneth == -1) { # if -1, solution of the ith column is 0
        soln = append(soln, 0)
      } else { # else, solution of the ith column is the value at the oneth row, last column
        soln = append(soln, tableau[oneth, ncol(tableau)])
      }
    }
  } else { # if it's a minimization problem, solution depends on the last row
    soln = tableau[nrow(tableau),] # solution is tableau's last row
    soln = soln[-(ncol(tableau)-1)] # remove second to the last value (z's value is the last)
  }
  
  basic.solution <- matrix( # store basic solution here
    data=soln, # data of basic.solution is equal to soln
    nrow=1, # single row; as much cols as tableau minus 1 (ans col)
    ncol=ncol(tableau)-1,
    dimnames = list("Solution", colnames(tableau)[1:(ncol(tableau)-1)]) # colnames is colnames of tableau minus Ans
  )
  
  if(problem) { # if problem is true, create new matrix
    xVals <- soln[9:(length(soln)-1)] # xVals refer to the 9th value of soln up to the penultimate one. in the given problem, this refers to x1 to x15
    shipping.num <- matrix( # matrix to put number of items to ship per plant to warehouse
      data=xVals, # xVals (x1-x15) is used as data
      nrow=3, # 3 rows (plants)
      ncol=5, # 5 columns (warehouses)
      byrow=TRUE, # add by row, indicate dimnames
      dimnames=list(c("Denver", "Phoenix", "Dallas"),c("Sacramento", "Salt Lake City", "Albuquerque", "Chicago", "New York City"))
    ) # form list with tableau, basic.solution, opt.val (z value, last value in basic.sol), and shipping.num
    list <- list(final.tableau = tableau, basic.solution = basic.solution, opt.val = basic.solution[nrow(basic.solution),ncol(basic.solution)], shipping.num=shipping.num)
  } else { # if problem is false, omit shipping num from return data
    list <- list(final.tableau = tableau, basic.solution = basic.solution, opt.val = basic.solution[nrow(basic.solution),ncol(basic.solution)])
  }
  return(list)
}
#TEST CASES IN HANDOUT
# max <- c(7,11,1,0,0,0,0,77,10,8,0,1,0,0,0,80,1,0,0,0,1,0,0,9,0,1,0,0,0,1,0,6,-150,-175,0,0,0,0,1,0)
# min <- c(1,2,1,0,0,4,7,6,0,1,0,20,-14,-20,0,0,1,0)
# init <- matrix(
#   data = max,
#   nrow = 5,
#   ncol = 8,
#   byrow = TRUE,
#   dimnames = list(c(1:5), c("x1","x2","S1","S2","S3","S4", "Z","Ans"))
# )
# print(simplex(init,TRUE,FALSE))
# 
# init <- matrix(
#   data = c(-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,-310,
#            0,0,0,0,0,-1,-1,-1,-1,-1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-260,
#            0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,0,0,1,0,0,0,0,0,0,-280,
#            1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,180,
#            0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,80,
#            0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,200,
#            0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,160,
#            0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,220,
#            -10,-8,-6,-5,-4,-6,-5,-4,-3,-6,-3,-4,-5,-5,-9,0,0,0,0,0,0,0,0,1,0),
#   nrow = 9,
#   ncol = 25,
#   byrow = TRUE,
#   dimnames = list(c(1:9), c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","S1","S2","S3","S4","S5","S6","S7","S8","Z","Ans"))
# )
# print(simplex(init,FALSE,TRUE))
# 
# init <- matrix(
#   data = c(-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,-200,
#            0,0,0,0,0,-1,-1,-1,-1,-1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-200,
#            0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,0,0,1,0,0,0,0,0,0,-200,
#            1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,100,
#            0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,100,
#            0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,100,
#            0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,100,
#            0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,100,
#            -5,-6,-7,-8,-9,-6,-7,-8,-9,-10,-3,-5,-7,-11,-13,0,0,0,0,0,0,0,0,1,0),
#   nrow = 9,
#   ncol = 25,
#   byrow = TRUE,
#   dimnames = list(c(1:9), c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","S1","S2","S3","S4","S5","S6","S7","S8","Z","Ans"))
# )
# print(simplex(init,FALSE,TRUE))
# 
# init <- matrix(
#   data = c(-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,-1400,
#            0,0,0,0,0,-1,-1,-1,-1,-1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-400,
#            0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,0,0,1,0,0,0,0,0,0,-200,
#            1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,431,
#            0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,332,
#            0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,350,
#            0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,450,
#            0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,400,
#            -30,-29,-31,-35,-33,-26,-24,-23,-25,-27,-11,-13,-15,-20,-17,0,0,0,0,0,0,0,0,1,0),
#   nrow = 9,
#   ncol = 25,
#   byrow = TRUE,
#   dimnames = list(c(1:9), c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","S1","S2","S3","S4","S5","S6","S7","S8","Z","Ans"))
# )
# print(simplex(init,FALSE,TRUE))
# 
# init <- matrix(
#   data = c(-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,-100,
#            0,0,0,0,0,-1,-1,-1,-1,-1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-100,
#            0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,0,0,1,0,0,0,0,0,0,-100,
#            1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,20,
#            0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,20,
#            0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,20,
#            0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,20,
#            0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,20,
#            -5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,0,0,0,0,0,0,0,0,1,0),
#   nrow = 9,
#   ncol = 25,
#   byrow = TRUE,
#   dimnames = list(c(1:9), c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","S1","S2","S3","S4","S5","S6","S7","S8","Z","Ans"))
# )
# print(simplex(init,FALSE,TRUE))
# 
# init <- matrix(
#   data = c(-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,-50,
#            0,0,0,0,0,-1,-1,-1,-1,-1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-50,
#            0,0,0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,0,0,1,0,0,0,0,0,0,-50,
#            1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,20,
#            0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,25,
#            0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,90,
#            0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,60,
#            0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,70,
#            -30,-29,-31,-35,-33,-26,-24,-23,-25,-27,-11,-13,-15,-20,-17,0,0,0,0,0,0,0,0,1,0),
#   nrow = 9,
#   ncol = 25,
#   byrow = TRUE,
#   dimnames = list(c(1:9), c("x1","x2","x3","x4","x5","x6","x7","x8","x9","x10","x11","x12","x13","x14","x15","S1","S2","S3","S4","S5","S6","S7","S8","Z","Ans"))
# )
# print(simplex(init,FALSE,TRUE))