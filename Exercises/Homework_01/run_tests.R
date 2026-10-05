tests_dir <- "~/Master CBB/Computational Biology/Exercises/Homework_01/student_test_suite" # replace with local path to the student_test_suite folder
if(tests_dir == "path/to/tests") stop("tests_dir needs to be set to a proper path")

library("RUnit")

testsuite <- defineTestSuite("HW", tests_dir)
currentdir <- getwd()
setwd(tests_dir)

out <- runTestSuite(testsuite, verbose=TRUE)
printTextProtocol(out)

setwd(currentdir)
