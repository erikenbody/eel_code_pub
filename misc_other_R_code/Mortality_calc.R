#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se


#Starting frequencies
p0 <- 0.5
q0 <- 0.5

#Final frequencies - if survival rates are 1 and d, respectively
#p1 <- p0/(p0 + q0*d)
#q1 <- q0*d/(p0 + q0*d)

#p1 - q1 = 0.1 desired DAF
#(1-d)/(1+d) = 0.1
#(1-d) = 0.1 + 0.1d
#10 - 10d = 1 + d
#11d = 9
d <- 9/11
p1 <- p0/(p0 + q0*d)
q1 <- q0*d/(p0 + q0*d)

pop0 <- p0 + q0
pop1 <- (p0 + q0*d)

#Assuming independent loci, each selection pressure should reduce the population by the same proportion, thus we get:
pop1_n10 <- pop1^10

#g1 <-0.5 - het freq
#g2 <- 0.25 - hom freq
#q0*d <- 0.5*g1 + g2*f2
#(q0*d - 0.25)/0.25 = f2
f2 <- (q0*d - 0.25)/0.25 

#fitness map
#0  = 1; 1 = 1;  2 = f2
fit_vec <- c(1,1,f2)

pop_df <- data.frame(id = 1:10000, stringsAsFactors = F)
for(i in 1:10){
  pop_df[,paste0("loc", i)] <- sample(rep(c(0,1,2), times = c(dim(pop_df)[1]/4, dim(pop_df)[1]/2, dim(pop_df)[1]/4)), size = dim(pop_df)[1])
  pop_df[pop_df[,paste0("loc", i)] == 0,paste0("loc", i,"_t1")] <- T
  pop_df[pop_df[,paste0("loc", i)] == 1,paste0("loc", i,"_t1")] <- T
  pop_df[pop_df[,paste0("loc", i)] == 2,paste0("loc", i,"_t1")] <- sample(c(T,F), prob = c(f2, 1-f2), replace = T, size = sum(pop_df$loc1 == 2))
}
t1_cols <- grep("_t1", names(pop_df))
survior_ids <- rowSums(pop_df[, t1_cols]) == length(t1_cols)
q1_vec <- numeric()
for(i in 1:10){
  tmp_table <- table(pop_df[survior_ids,t1_cols[i]-1])
  q1_vec[i] <- (tmp_table[2] * 0.5 + tmp_table[3])/sum(tmp_table)
}


