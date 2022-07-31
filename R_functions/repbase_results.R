repbase_results <- function(repbase.csv){
  csv <- read.csv(repbase.csv, header = T, sep = ",")
  with_max_score <- csv %>% group_by(Name_seq) %>% mutate(max = max(Score))
  with_is_max <- with_max_score %>% group_by(Name_seq) %>% mutate(is_max = (Score == max))
  count <- with_is_max %>%  group_by(Name_seq) %>% count(Name_seq)
  
  df <- inner_join(with_is_max, count)
  df <- df %>% group_by(Name_seq) %>% mutate(non_max_sum = sum(Score[is_max == FALSE]))
  df <- df %>% group_by(Name_seq) %>% mutate(non_max_avg = non_max_sum/(n-1))
  
  max_classes <- df %>% filter(is_max == TRUE)
  #adding in the threshold
  df <- max_classes %>% mutate(by_2.5 = (Score >= 2.5 * non_max_avg))
  single_copy <- df %>% filter(n == 1)
  repbase_results <- df %>% filter(by_2.5 == TRUE)
  
  res <- rbind.data.frame(single_copy, repbase_results)
  
  results <- data.frame("Name" = res$Name_seq, 
                        "Length" = res$To_seq - res$From_seq,
                        "Class" = res$Class,
                        "Sim" = res$Sim,
                        "Score" = res$Score,
                        "non_max_avg" = res$non_max_avg,
                        "unique_hits" = res$n,
                        "by>2.5" = res$by_2.5)
  results
}
