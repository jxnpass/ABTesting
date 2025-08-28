# Huge credits to @BillVenables
# "https://rdrr.io/github/BillVenables/WWRCourse/src/R/roundRobin.R"
# as this code comes straight from their work.

#### function to schedule a round robin tournament

#' Generate a round robin tournament schedule
#'
#' Given a list of n team names, generate a program of n-1 Rounds
#' where each team plays each other team precisely once.  If n is
#' odd the teams are augmented by a dummy team <Bye>.
#'
#' @param teams Either an interger specifing the number of teams,
#'        or a character string vector giving their names
#' @param alphabetical logical: should the teams be alphabetically
#'        ordered, if necessary?
#' @param reorder logical: within each round should the games be
#'        listed in alphabetical order of the "Home" team?
#'
#' @return A 3-dimensional (n/2 x 2 x (n-1)) array giving the
#'        entire tournament
#' @export
#'
#' @examples
#' (Season2019 <- round_robin(NRL))
#' summary(Season2019)
#' summary(Season2019, "travel")
round_robin <- function(teams, alphabetical = TRUE, reorder = FALSE) {
  # If a single numeric value for teams, generate default team names
  if(is.numeric(teams) && length(teams) == 1 && teams > 0) {
    teams <- paste0("Team", 1:ceiling(teams))
  }
  
  stopifnot(is.character(teams) && length(teams) > 0)
  
  # Check for duplicates in team names
  if(any(duplicated(teams))) 
    stop("duplicated team names are not allowed")
  
  # Sort teams alphabetically if specified
  if(alphabetical) teams <- sort(teams)
  
  # Add a "Bye" if the number of teams is odd
  odd <- length(teams) %% 2 == 1
  if(odd) 
    teams <- c("<Bye>", teams)
  
  n <- length(teams)
  
  # Generate the initial round by pairing teams
  round <- cbind(teams[1:(n/2)], teams[n:(n/2+1)])
  colnames(round) <- c("Home", "Away")
  rownames(round) <- paste("Match", 1:nrow(round))  # Simplified match labeling
  
  # Initialize result list with the first round
  Res <- structure(vector("list", n-1),
                   names = paste("Round", seq_len(n-1)))  # Simplified round naming
  Res[[1]] <- round
  
  # Generate remaining rounds
  if(n > 2) {
    ij <- rbind(cbind(2:nrow(round), 1), cbind(nrow(round):1, 2))
    ji <- rbind(ij[-1,,drop = FALSE], ij[1,,drop = FALSE])
    
    for(i in 2:(n-1)) {
      round[ji] <- round[ij]
      r <- round
      if(!odd && (i %% 2 == 0)) 
        r[1, 1:2] <- r[1, 2:1]
      Res[[i]] <- r
    }
  }
  
  # Reorder the rounds if specified
  if(reorder) {
    for(j in seq_along(Res)) 
      Res[[j]][] <- as.vector(Res[[j]][order(Res[[j]][, 1]), ])
  }
  
  # Assign class to the result
  class(Res) <- "round_robin"
  
  return(Res)
}


#' @rdname round_robin
#' @export
summary.round_robin <- function(object, type = c("venue", "travel"), ...) {
  type <- match.arg(type)
  if(type == "travel") {
    rounds <- format(seq_along(object))
    for(i in seq_along(object)) {
      object[[i]] <- cbind(Round = rounds[i], object[[i]])
    }
  }
  g <- do.call(rbind, object)
  if(any(byes <- grepl("^<Bye> *$", g))) {
    dim(byes) <- dim(g)
    i <- which(byes, arr.ind = TRUE)[, "row"]
    g <- g[-i,]
  }
  if(type == "venue") {
    f <- factor(col(g), levels = 1:2, labels = colnames(g))
    table(Team = g, Venue = f)
  } else {
    # g <- data.frame(g, stringsAsFactors = FALSE) %>%
    #   gather(key = Venue, value = Team, Home, Away) %>%
    #   within(Venue <- format(substring(Venue, 0, 1),
    #                          justify = "right",
    #                          width = nchar(rounds[1]))) %>%
    #   spread(key = Round, value = Venue)
    Home <- cbind(g[, c("Round", "Home")], Venue = "Home")
    Away <- cbind(g[, c("Round", "Away")], Venue = "Away")
    g <- rbind(Home, Away)
    colnames(g)[2] <- "Team"
    g <- within(data.frame(g, stringsAsFactors = FALSE), {
      Venue <- format(substring(Venue, 0, 1),
                      justify = "right",
                      width = nchar(rounds[1]))
    })
    teams <- with(g, sort(unique(Team)))
    rounds <- with(g, sort(unique(Round)))
    out <- matrix(NA_character_, length(teams), length(rounds))
    dimnames(out) <- list(Team = teams, Round = rounds)
    ij <- with(g, cbind(match(Team, teams), match(Round, rounds)))
    out[ij] <- g$Venue
    # out <- as.matrix(g[, -1])
    out[is.na(out)] <- format("*", width = nchar(rounds[1]),
                              justify = "right")
    # rownames(out) <- g[["Team"]]
    noquote(out)
  }
}

to_df <- function(rr) {
  
  rr_df <- tibble()
  
  for (round_name in names(rr)) {
    round_data <- as_tibble(rr[[round_name]])
    round_data$Round <- rep(round_name, each = nrow(round_data))
    round_data$Match <- paste("Match", 1:nrow(round_data))
    rr_df <- bind_rows(rr_df, round_data)
  }
  
  return(rr_df)
}

stack_matchups <- function(matchups) {
  matchups %>%
    imap_dfr(~ .x %>% mutate(Pool = .y)) %>% 
    mutate(ID = paste(Match, "-", Pool))
}

to_sheets <- function(matchups, link, stack){
  
  if (stack == "BOTH") {
    for (pool in pools){
      matchups[[pool]]$Pool <- pool
      matchups[[pool]]$ID <- paste(matchups[[pool]]$Match, "-", pool)
      sn <- paste("Pool", pool, "Schedule")
      sheet_write(data = matchups[[pool]], ss = link, sheet = sn)
    }
    stacked_matchups <- stack_matchups(matchups)
    sn <- "Pool Full Schedule"
    sheet_write(data = stacked_matchups, ss = link, sheet = sn)
  }
  
  else if (!stack) {
    for (pool in pools){
      matchups[[pool]]$Pool <- pool
      matchups[[pool]]$ID <- paste(matchups[[pool]]$Match, "-", pool)
      sn <- paste("Pool", pool, "Schedule")
      sheet_write(data = matchups[[pool]], ss = link, sheet = sn)
    }
  }
  
  else if (stack) {
    stacked_matchups <- stack_matchups(matchups)
    sn <- "Pool Full Schedule"
    sheet_write(data = stacked_matchups, ss = link, sheet = sn)
  }
  
  print("Done")
} 

