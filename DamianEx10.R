source("DamianEx08.R") # import for qsi and simplex
source("DamianEx09.R")

library(shiny) # import for gui
library(shinythemes)

ui <- fluidPage(theme = shinytheme("cerulean"),
                tags$head( # styling for prompt; add border, padding color, etc.
                  tags$style(HTML("
                  #prompt {
                    padding: 5px;
                    border:1px solid;
                    border-radius: 5px;
                    border-color: #bfbfbf;
                    background-color: white;
                    color:black;
                  }"))),
                navbarPage( # add ribbons at the top of page
                  "CMSC150 Integration - Damian", # title
                  tabPanel("Quadratic Spline Interpolation", # first tab is for qsi
                           sidebarPanel( # at the side are the input spaces
                             textInput("x", "Enter x values:", ""),
                             textInput("y", "Enter y values:", ""),
                             numericInput("z", "Enter x value to approximate:", ""),
                           ),
                           mainPanel( # in the center is the output
                             verbatimTextOutput("qsi.fxns", placeholder=TRUE), # automatically prints output based on x,y,z above
                             
                           ) 
                           
                  ),
                  tabPanel("Simplex Programming", # second tab is for simplex
                            sidebarPanel( # at the side are the input spaces
                              fluidRow( # there are two select inputs beside each other
                                column(6, selectInput("isMax", label = NULL, width='100%', # the first is for max/min; default is max
                                            choices = list("Maximize" = TRUE, "Minimize" = FALSE), 
                                            selected = TRUE)),
                                column(6, selectInput("problem", label = NULL, width='100%', # the second is for problem/new; default is new
                                            choices = list("Given Problem" = TRUE, "New Problem" = FALSE), 
                                            selected = FALSE)),
                              ),
                              textOutput("prompt"), # print prompt based on isMax
                              br(), # linebreak, then ask for tabluea vals. print placeholder to give example
                              textAreaInput("tableauVals", label="Enter  constraints and the objective function:", resize='both', placeholder="1,2,4,\n7,6,20,\n-14,-20,0"),
                            ),
                            mainPanel( # at the center is the automatic output
                              verbatimTextOutput("result", placeholder=TRUE),
                            ) 
                )
           )
)

server <- function(input, output) {
  
  output$qsi.fxns <- renderPrint({ # prints output for qsi
    tryCatch( # catch input errors
      expr = {
        x=eval(parse(text=paste("c(",input$x,")",sep=""))) # enclose input x in c(), parse, then eval to make it a vector
        y=eval(parse(text=paste("c(",input$y,")",sep="")) )# enclose input y in c(), parse, then eval to make it a vector
        data=list(xVals=x,yVals=y) # combine x and y in a list
        result = poly.qsi(data, input$z) # solve qsi with data and z as inputs
        print("Quadratic Polynomials:")
        for(i in result$qsi.fxns) { # print each polynomials
          print(paste(deparse(i)[1], deparse(i)[2], sep="")) # paste to combine "function (x)" with "x1 + ..."
        }
        print("-----------------------------------")
        print(paste("Approximate y-value at ", input$z, ": ", result$y, sep="")) # print estimated y value
        },
      error = function(e){ # if an error is encountered, print a prompt
        return("Enter valid inputs")
      }
    )
    })
  
  output$prompt <- renderText({ # prints a prompt depending on input isMax
    if (input$isMax) {
      return("Maximization: All constraints should have a <= inequality. If there are >= constraints, multiply the values by -1 to flip the inequality. Enter only the coeficients (eg. x1 + 8x3 <= 9 → 1, 0, 8, 9). For the objective function, the coefficients are multiplied by -1 while the RHS is 0.")
    } else {
      return("Minimization: All constraints should have a >= inequality. If there are <= constraints, multiply the values by -1 to flip the inequality. Enter only the coeficients (eg. x1 + 8x3 <= 9 → 1, 0, 8, 9). For the objective function, the coefficients are multiplied by -1 while the RHS is 0.")
    }
  })
  
  output$result <- renderPrint({ # prints output for simplex
    tryCatch( # catch input errors
      expr = {
        funcs = unlist(strsplit(input$tableauVals, "\n")) # make a vector composed of every input line (whole input was split by \n)
        numConstraints = length(funcs)-1 # number of constraints depends on number of funcs - 1 (because of the objective func)
        funcs[numConstraints+1] = paste(funcs[numConstraints+1], ",", sep="") # add a comma at the end of the last line for formatting
        numVars = length(unlist(strsplit(funcs[1], ",")))-1 # number of variables depend on number of values when a line is split by a comma - 1 (the rhs)

        tableauVals = c() # store all values combined (plus slacks) here
        slacks = numeric(length=numConstraints+1) # per line, there would be slacks. Initially all 0, length depends on number of constraints + 1 (for Z)
        for (i in 1:length(funcs)) { # loop through every func
          coeffs = unlist(strsplit(funcs[i], ",")) # store coefficients here (split by comma)
          xVars = "" # store x vars here (exclude last value in line since that is for rhs)
          for(j in 1:numVars) { # append "x, " in xVars for every x
            xVars=paste(xVars, coeffs[j], ", ", sep="")
          }

          slacks[i] = 1 # the ith value in slacks must be 1
          slacksStr = "" # store slacks for the ith func here
          for(j in 1:length(slacks)) { # loop through every slack, append slack's jth value and a comma to slacksStr
            slacksStr = paste(slacksStr, slacks[j], ", ", sep="")
          }
          # the iRow is composed of xVars, followed by slacks, then the last value of coeffs (rhs). enclose in c(), parse, then eval to make it a vector
          iRow = eval(parse(text=paste("c(",  xVars, slacksStr, coeffs[length(coeffs)] ,")",sep=""))) 
          tableauVals = append(tableauVals, iRow) # append iRow to tableauVals
          slacks = numeric(length=numConstraints+1) # reset slacks to all 0 for the next func
        }
        
        columnNames <- c() # store column name of matrix to set up here
        for(i in 1:numVars) {columnNames = append(columnNames, paste("x", i, sep=""))} # loop through vars, append Xi
        for(i in 1:numConstraints) {columnNames = append(columnNames, paste("S", i, sep=""))} # loop through slacks, append Si
        columnNames = append(columnNames, c("Z", "Ans")) # append Z and Ans

        init <- matrix( # set up matrix here
          data = tableauVals, # data corresponds to tableauVals
          nrow = numConstraints+1, # number of rows depends on number of constraints + 1 (for objective function)
          ncol = numVars + numConstraints + 2, # number of column depends on numVars + numConstraints + 2 (for Z and Ans)
          byrow = TRUE, # add by row
          dimnames = list(c(1:(numConstraints+1)), columnNames) # rownames is just a numbering, colnames is from columnNames
        )
        return(simplex(init, eval(parse(text=input$isMax)), eval(parse(text=input$problem)))) # solve thru simplex, use init, input$isMax, and input$problem
      },
      error = function(e) { # if an error is encountered, print a prompt
        return("Enter valid inputs")
      }
    )
  })
}

shinyApp(ui = ui, server = server) # combine ui and server