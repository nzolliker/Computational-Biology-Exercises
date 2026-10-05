#####################################################
#######        COMPUTATIONAL BIOLOGY         ########
#######             HOMEWORK 1               ########
#####################################################
#                                                   #
# Implement the pairwise alignment algorithms       #
# Needleman-Wunsch and Smith-Waterman.              #
#                                                   #
#####################################################
#####################################################

# In all functions the following parameters are the same:
# seqA: the first sequence to align
# seqB: the second sequence to align
# score_gap: score for a gap
# score_match: score for a character match
# score_mismatch: score for a character mismatch
# local: (logical) True if alignment is local, False otherwise

init_score_matrix = function(nrow, ncol, local, score_gap) {
    # Initialize the score matrix with zeros.
    # If the alignment is global, the leftmost column and the top row will have incremental gap scores,
    # i.e. if the gap score is -2 and the number of columns is 4, the top row will be [0, -2, -4, -6].
    # nrow: (numeric) number of rows in the matrix
    # ncol: (numeric) number of columns in the matrix

    # ???
    score_matrix <- matrix(0, nrow, ncol)
    if (local == FALSE){
      for (j in seq(2, ncol, 1)){
        score_matrix[1,j] <- (j-1)*score_gap
      }
      for (i in seq(2, nrow, 1)){
        score_matrix[i,1] <- (i-1)*score_gap
      }
    }
    # Return the initialized empty score matrix
    # score_matrix: (numeric) nrow by ncol matrix
    return(score_matrix)
}

init_path_matrix = function(nrow, ncol, local) {
    # Initialize the path matrix with empty values ("").
    # Additionally, for GLOBAL alignment (i.e. local==FALSE), make the first row
    # have "left" on all positions except 1st, and make the first column
    # have "up" on all positions except 1st.
    # nrow: (numeric) number of rows in the matrix
    # ncol: (numeric) number of columns in the matrix

    # ???
  path_matrix <- matrix("", nrow, ncol)
  if (local == FALSE){
    for (j in seq(2, ncol, 1)){
      path_matrix[1,j] <- "left"
    }
    for (i in seq(2, nrow, 1)){
      path_matrix[i,1] <- "up"
    }
  }

    # Return the initialized empty path matrix
    # path_matrix: (character) nrow by ncol matrix
    return(path_matrix)
}

get_best_score_and_path = function(row, col, nucA, nucB, score_matrix, score_gap, score_match, score_mismatch, local) {
    # Compute the score and the best path for a particular position in the score matrix
    # nucA: (character) nucleotide in sequence A
    # nucB: (character) nucleotide in sequence B
    # row: (numeric) row-wise position in the matrix
    # col: (numeric) column-wise position in the matrix
    # score_matrix: (double) the score_matrix that is being filled out

    # ???
  if (nucA == nucB) {
    case_1 <- score_matrix[row-1,col-1] + score_match
  } else {
    case_1 <- score_matrix[row-1,col-1] + score_mismatch
  }
  case_2 <- score_matrix[row, col-1] + score_gap
  case_3 <- score_matrix[row-1, col] + score_gap 
  
  if (local == TRUE){
    paths <- c("diag", "left", "up", "-")
    cases <- c(case_1, case_2, case_3, 0)
    for (i in seq(1,4,1)){
      if (cases[i] == max(cases)){
        score <- cases[i]
        path <- paths[i]
      }
    }
  } else {
    paths <- c("diag", "left", "up")
    cases <- c(case_1, case_2, case_3)
    for (i in seq(1,3,1)){
      if (cases[i] == max(cases)){
        score <- cases[i]
        path <- paths[i]
      }
    }
  }
  
  
    # Return the best score for the particular position in the score matrix
    # In the case that there are several equally good paths available, return any one of them.
    # score: (numeric) best score at this position
    # path: (character) path corresponding to the best score, one of ["diag", "up", "left"] in the global case and of ["diag", "up", "left", "-"] in the local case
    return(list("score"=score, "path"=path))
}

fill_matrices = function(seqA, seqB, score_gap, score_match, score_mismatch, local, score_matrix, path_matrix) {
    # Compute the full score and path matrices
    # score_matrix: (numeric)  initial matrix of the scores
    # path_matrix: (character) initial matrix of paths

    # ???
  for (i in seq(2, nchar(seqA)+1, 1)){
    for (j in seq(2, nchar(seqB)+1, 1)){
      nucA <- substr(seqA, i-1, i-1)
      nucB <- substr(seqB, j-1, j-1)
      row <- i
      col <- j
      x <- get_best_score_and_path(row, col, nucA, nucB, score_matrix, score_gap, score_match, score_mismatch, local) 
      score_matrix[i,j] <- x$score
      path_matrix[i,j] <- x$path
    }
  }

    # Return the full score and path matrices
    # score_matrix: (numeric) filled up matrix of the scores
    # path_matrix: (character) filled up matrix of paths
    return(list("score_matrix"=score_matrix, "path_matrix"=path_matrix))
}

get_best_move = function(nucA, nucB, path, row, col) {
    # Compute the aligned characters at the given position in the score matrix and return the new position,
    # i.e. if the path is diagonal both the characters in seqA and seqB should be added,
    # if the path is up or left, there is a gap in one of the sequences.
    # nucA: (character) nucleotide in sequence A
    # nucB: (character) nucleotide in sequence B
    # path: (character) best path pre-computed for the given position
    # row: (numeric) row-wise position in the matrix
    # col: (numeric) column-wise position in the matrix

    # ???
  if (path == "diag"){
    char1 <- nucA
    char2 <- nucB
    newrow <- (row-1)
    newcol <- (col-1)
  }
  if (path == "up"){
    char1 <- nucA
    char2 <- "-"
    newrow <- (row-1)
    newcol <- col
  }
  if (path == "left"){
    char1 <- "-"
    char2 <- nucB
    newrow <- row
    newcol <- (col-1)
  }
  if (path == "-"){
    char1 <- ""
    char2 <- ""
    newrow <- 0
    newcol <- 0
  }

    # Return the new row and column and the aligned characters
    # newrow: (numeric) row if gap in seqA, row - 1 otherwise
    # newcol: (numeric) col if gap in seqB, col - 1 otherwise = up
    # char1: (character) '-' if gap in seqA, appropriate character if a match
    # char2: (character) '-' if gap in seqB, appropriate character if a match = up
    return(list("newrow"=newrow, "newcol"=newcol, "char1"=char1, "char2"=char2))
}

get_best_alignment = function(seqA, seqB, score_matrix, path_matrix, local) {
    # Return the best alignment from the pre-computed score matrix
    # score_matrix: (numeric) filled up matrix of the scores
    # path_matrix: (character) filled up matrix of paths

    # ???
  
  alignment <- c("", "")
  if (local) {
    score <- max(score_matrix)
    ali <- which(score_matrix == score, arr.ind = TRUE)[1,]
    row <- as.integer(ali[1])
    col <- as.integer(ali[2])
    while (score_matrix[row, col] > 0){
      nucA <- substr(seqA, row-1, row-1)
      nucB <- substr(seqB, col-1, col-1)
      path <- path_matrix[row, col]
      values <- get_best_move(nucA, nucB, path, row, col)
      row <- values$newrow
      col <- values$newcol
      alignment[1] <- paste0(values$char1, alignment[1])
      alignment[2] <- paste0(values$char2, alignment[2])      
    }
  } else {
    row <- nchar(seqA)+1
    col <- nchar(seqB)+1
    score <- score_matrix[row, col]
    while (row > 1 || col > 1) {
      nucA <- substr(seqA, row-1, row-1)
      nucB <- substr(seqB, col-1, col-1)
      path <- path_matrix[row, col]
      values <- get_best_move(nucA, nucB, path, row, col)
      row <- values$newrow
      col <- values$newcol
      alignment[1] <- paste0(values$char1, alignment[1])
      alignment[2] <- paste0(values$char2, alignment[2])
    }
  }
    # Return the best score and alignment (or one thereof if there are multiple with equal score)
    # score: (numeric) score of the best alignment
    # alignment: (character) the actual alignment in the form of a vector of two strings
    return(list("score"=score, "alignment"=alignment))
}

align = function(seqA, seqB, score_gap, score_match, score_mismatch, local) {
    # Align the two sequences given the scoring scheme
    # For testing purposes, use seqA for the rows and seqB for the columns of the matrices
  
    # Initialize score and path matrices
    # ???
    nrow <- nchar(seqA)+1
    ncol <- nchar(seqB)+1
    score_matrix <- init_score_matrix(nrow, ncol, local, score_gap)
    path_matrix <- init_path_matrix(nrow, ncol, local)
  
    # Fill in the matrices with scores and paths using dynamic programming
    # ???
  y <- fill_matrices(seqA, seqB, score_gap, score_match, score_mismatch, local, score_matrix, path_matrix)
  
    # Get the best score and alignment (or one thereof if there are multiple with equal score)
    # ???
  result <- get_best_alignment(seqA, seqB, y$score_matrix, y$path_matrix, local)
    
    # Return the best score and alignment (or one thereof if there are multiple with equal score)
    # Returns the same value types as get_best_alignment
    return(result)
}

test_align = function() {
    seqA = "TCACACTAC"
    seqB = "AGCACAC"
    score_gap = -2
    score_match = +3
    score_mismatch = -1
    local = F
    result = align(seqA, seqB, score_gap, score_match, score_mismatch, local)
    print(result$alignment)
    print(result$score)
}

test_align()

