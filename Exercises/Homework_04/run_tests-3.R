tests_dir <- "~/Master CBB/Computational Biology/Exercises/Homework_04/student_test_suite" # replace with local path to the student_test_suite folder

if(!dir.exists(tests_dir)) {
  stop("Please set \"tests_dir\" as the local path to your student_test_suite folder.")
}

library("RUnit")

original_dir <- getwd()
setwd(tests_dir)

testsuite <- defineTestSuite("HW", ".", rngKind=RNGkind()[1], rngNormalKind=RNGkind()[2])
out <- runTestSuite(testsuite)
printTextProtocol(out)

setwd(original_dir)
