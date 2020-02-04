#Ex 2.1
sample(c("ii", "iB", "iA", "BB", "AA", "AB"), replace = TRUE, size = 1)

#Ex 2.2
offspringG <- function(Parent1, Parent2){
  Allele1 <- sample(unlist(strsplit(Parent1, split = "")), replace = TRUE, size=1)
  Allele2 <- sample(unlist(strsplit(Parent2, split = "")), replace = TRUE, size=1)
  return(paste(Allele1, Allele2, sep = ""))
}
offspringG("AA", "AB")
offspringP <- function(Genotype){
  Phenotype <- ""
  if (Genotype == "ii"){
    Phenotype <- "Type O"
    }
  if (Genotype == "AB" | Genotype == "BA"){
    Phenotype <- "Type AB"
    }
  if (Genotype == "iB" | Genotype == "Bi" | Genotype == "BB"){
    Phenotype <- "Type B"
  }
  if (Genotype == "iA" | Genotype == "Ai" | Genotype == "AA"){
    Phenotype <- "Type A"
  }
  return(Phenotype)
}

#Ex2.3
offspingiAiA <- c()
for(i in 1:1000){
  offspingiAiA[i] <- offspringG("iA", "iA")
}
offspingiAiA == "ii"

sum(offspingiAiA == "ii")/length(offspingiAiA)
sum(offspingiAiA != "AA")/length(offspingiAiA)

#Ex2.4
secondGen <- c()
for(i in 1:1000){
  A <- offspringG("iA", "iB")
  B <- offspringG("AB", "AB")
  C <- offspringP(offspringG(A, B))
  secondGen[i] <- C
}
sum(secondGen == "Type AB")/length(secondGen)


#Ex4.1
ACGsim <- function(reps, length){
  count <- 0
  for(i in 1:reps){
    seq <- sample(c("A", "C", "T", "G"), replace =TRUE, size = length)
    for(i in 3:length(seq)){
      if(seq[i-2] == "A" & seq[i-1] == "C" & seq[i] == "G"){
        count <- count + 1
      }
    }
  }
  return(count/reps)
}

ACGsim(10000, 5)
system.time(ACGsim(10000, 5))
x <- system.time(JonathanBaptism(10000, 5))
JonathanBaptism(10000, 5)

JonathanBaptism <- function(reps, length = 10){
  count <- 0
  for (i in 1:reps){
    seq <- paste(sample(c("A", "C", "T", "G"), replace =TRUE, size = length), collapse = "")
    if(grepl("ACG", seq, fixed=TRUE)){
      count <- count + 1
    }
  }
  return(count/reps)
}

#Ex 4.2




Gsim <- function(reps, length = 10){
  miss <- 0
  for(i in 1:reps){
    seq <- sample(c("A", "C", "T", "G"), replace =TRUE, size = length)
    if ("G" %in% seq){
      miss <- miss + 1
    }
  }
  count <- reps - miss
  return(count/reps)
}
Gsim(10000)

#Ex 4.3
simG <- c()
for(i in 1:20){
  simG[i] <- Gsim(10000, i)
}
pG <- c()
for(i in 1:20){
  pG[i] <- 0.75**i
}
plot(1:20, simG, col = "blue", pch = 0)
points(1:20, pG, col = "red", pch = 1)
legend("topright",legend = c("simulation", "probability"), col = c("blue", "red"), pch = c(0, 1))

#Ex 5.1
isidata <- read.table("assignment1/neuronspikes.txt", col.names = "isi")
hist(isidata$isi, breaks = 50, probability = TRUE)

#Ex 5.2
for(i in seq(0.5, 1.5, 0.1)){
  curve(dexp((x), rate = i), from = 0, to = 5, col= (i*10)-3, add = TRUE)
}
legend("topright", legend=c(seq(0.5,1.5,0.1)), col=c(2:length(seq(0.5,1.5,0.1))), lty = 1)
#Which is best fit? Purple maybe? not sure

#Ex 5.3

#Ex 6.1
cells <- read.csv("assignment1/cell_types.csv", na.strings = "") #Lots of missing values, na.strings turns them into NA values

#Ex.6. 2
#The donor species column indicates which values are humans and which are mice
sum(cells$donor__species == "Homo Sapiens")/length(cells$donor__species)

#Ex6.3
hist(cells$ef__peak_t_ramp[cells$donor__species == "Homo Sapiens"], breaks = 50, probability = TRUE)
hist(cells$ef__peak_t_ramp[cells$donor__species == "Mus musculus"], breaks = 50, probability = TRUE)

#Ex 6.4
first_plot = FALSE
for(i in 1:3) {
  curve(
    dlnorm((x), sdlog = 0.6, meanlog = i),
    from = 0,
    to = 30,
    col = i + 1,
    add = first_plot
  )
  first_plot = TRUE
}
legend("topright", legend=c(1:3), col=c(2:4), lty = 1)

#Ex6.5
#2 is best fit

#Ex 6.6
hist(cells$ef__peak_t_ramp[cells$donor__species == "Homo Sapiens"], breaks = 50, probability = TRUE)
curve(dlnorm((x), sdlog = 0.6, meanlog = 2), from = 0, to = 30, col=3, add = TRUE)
#Fits alright I guess

#Ex 6.7
sum(cells$donor__species == "Homo Sapiens" & cells$donor__sex == "Male")
sum(cells$donor__species == "Homo Sapiens" & cells$donor__sex == "Female")

