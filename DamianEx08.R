source("DamianEx04.R") # import for gauss jordan

poly.qsi <- function(data, x) { # data contains x and y vectors while x is the point to estimate data in
  # sort x values ascendingly and match y values in data (through selection sort)
  for (i in 1:length(data$xVals)) { # outer loop is for every x value
    least = i # will store index of smallest x value among unsorted values
    for (j in i:length(data$xVals)) { # loop from i to end (start at i since values before it are sorted already)
      if (data$xVals[j] < data$xVals[least]) { # if jth value is less than value in index least
        least = j # set j as the new least
      }
    }
    temp = data$xVals[i] # swap ith value with the least one for x
    data$xVals[i] = data$xVals[least]
    data$xVals[least] = temp

    temp2 = data$yVals[i] # y is not sorted but must be swapped as well to maintain pairing with their x values
    data$yVals[i] = data$yVals[least]
    data$yVals[least] = temp2
  }
  
  n = length(data$xVals)-1 # number of intervals is number of data points - 1

  sys <- matrix( # matrix to put xi values in
    data = 0, # initially all 0
    nrow = (3*n), # 3n rows/equations
    ncol = (3*n)+1, # +1 for rhs // order is a1, a2,.. an, b1, b2, .. cn, rhs
    byrow = TRUE
  )
  
  currRow = 1 # keep count on what eq/row we are in
  for(i in 2:n) { # 1st Condition concerning internal knots
    sys[currRow, i-1] = data$xVals[i] * data$xVals[i] # ai-1 = xi^2
    sys[currRow, i-1+n] = data$xVals[i] # bi-1 = xi (add n to get right index)
    sys[currRow, i-1+(2*n)] = 1 # ci-1 = 1 (add 2n to get right index)
    sys[currRow, (3*n)+1] = data$yVals[i] # rhs (last column) = yi
    currRow = currRow + 1 # move to next row
    
    sys[currRow, i] = data$xVals[i] * data$xVals[i] # similar to above but w/out the -1
    sys[currRow, i+n] = data$xVals[i]
    sys[currRow, i+(2*n)] = 1
    sys[currRow, (3*n)+1] = data$yVals[i]
    currRow = currRow + 1
  }
  # 2nd Condition concerning end points
  # first point (dont bother with a1)
  sys[currRow, 1+n] = data$xVals[1] # b1 = x1
  sys[currRow, 1+(2*n)] = 1 # c1 = 1
  sys[currRow, (3*n)+1] = data$yVals[1] # rhs = y1
  currRow = currRow + 1 
  # end point
  sys[currRow, n] = data$xVals[n+1] * data$xVals[n+1] # an = (xn+1)^2
  sys[currRow, n+n] = data$xVals[n+1] # bn = xn+1
  sys[currRow, n+(2*n)] = 1 # cn = 1
  sys[currRow, (3*n)+1] = data$yVals[n+1] # rhs = yn+1
  currRow = currRow + 1
  
  for (i in 2:n) { # 3rd condition concerning first derivative of interior knots
    sys[currRow, i-1] = data$xVals[i]*2 # ai-1 = 2xi
    sys[currRow, i] = data$xVals[i]*-2 # ai = -2xi
    sys[currRow, i-1+n] = 1 # bi-1 = 1
    sys[currRow, i+n] = -1 # bi = -1
    currRow = currRow + 1
  }
  # 4th condition: a1 = 0
  sys <- sys[-(3*n),-1] # remove 1st column and last row since a1 = 0
  input <- list(augcoeffmatrix = sys, variables = 1:((3*n)-1)) # form input to gauss jordan function
  ans = GaussJordanElimination(input)$solutionSet # call function to solve matrix
  
  qsi.fxns = list() # store quadratic functions here
  for(i in 1:n) { # n intervals and functions
    str = paste("function (x) ", sep="") # make function in form of a string first
    if (i==1) { # if i=1, then a = 0 so no x^2 term
      str = paste(str, round(ans[n],4),"*x + ", round(ans[(2*n)],4),sep="") # append b*x + c to str
    } else { # else, append a*x^2 + b*x + c
      str = paste(str, round(ans[i-1],4), "*x^2 + ", round(ans[i-1+n],4), "*x + ", round(ans[i-1+(2*n)],4), sep="")
    }
    str = eval(parse(text=str)) # parse then eval to make it a function
    qsi.fxns = append(qsi.fxns, str) # append function to list
  }
  # find estimated value at x
  y = NA # set y as NA initially (if x is out of range, y will remain NA)
  for (i in 1:n) { # loop through every interval
    if (data$xVals[i] <= x && data$xVals[i+1] >= x) { # if x is >= xi and <= xi+1, we are at the right interval
      y = round(qsi.fxns[[i]](x),4) # pass x (user input) to the ith function
      break 
    }
  }
  
  result = list(qsi.fxns=qsi.fxns, y=y) # labeled list with qsi.fxns and y
  return(result)
  
}

# data = list(xVals=c(7, 9, 3, 4.5), yVals=c(2.5, 0.5, 2.5, 1))
# print(poly.qsi(data, 5))