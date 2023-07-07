
Watson = as.matrix(d[,6:9])
Crick  = as.matrix(d[,11:14])[,c(2,1,4,3)]

# take sqrt transform of the counts
# considering the counts are actually positively correlated

t = 60
k = t - sqrt(t)

Watson[Watson > t] = round(sqrt(Watson[Watson > t]) + k)
Crick[Crick > t] = round(sqrt(Crick[Crick > t]) + k)

# mutation rate
pm = 1/1000/3

# error rate
# in oocyte samples, error rate is set triple
# pe = 1/100/3
pe = 3/100/3

# total mis rate
p = pm + pe

# methylation rate/proportion 
pr.cg = 0.6        # CG content
pr.ncg = 1/100     # non-CG content


# transition prob of haploidy

PA = function(pr){c(1-3*pm-3*pe, 2*pm-pm*pr+pe, pm*pr+pe, pm+pe)}
PT = function(pr){c(pm+pe, 1-2*pm-pm*pr-3*pe, pm*pr+pe, pm+pe)}
PC = function(pr){c(pm+pe, pm+pe+(1-3*pm-3*pe)*(1-pr), (1-3*pm-3*pe)*pr, pm+pe)}
PG = function(pr){c(pm+pe, 2*pm-pm*pr+pe, pm*pr+pe, 1-3*pm-3*pe)}


# STATUS
HOMO = c('A', 'T', 'C', 'G')
HETER = c('AC', 'AG', 'AT', 'CG', 'CT', 'GT')
STATUS = c(HOMO, HETER)

# prior

ps = c((1-3*p)^2, p^2, 2*p*(1-3*p))

pri.A = ps[c(1,2,2,2,3,3,3,2,2,2)]
pri.T = ps[c(2,1,2,2,2,2,3,2,3,3)]
pri.C = ps[c(2,2,1,2,3,2,2,3,3,2)]
pri.G = ps[c(2,2,2,1,2,3,2,3,2,3)]

pris = list(pri.A, pri.T, pri.C, pri.G)


postp <- function(ref = 'A', cg = TRUE, Watson = 1:4, Crick = 1:4) {
  
  if(ref == 'N') return(NA)
  
  if(cg) pr = pr.cg
  else pr = pr.ncg
  
  # prior
  
  theta = pris[[which.max(c('A', 'T', 'C', 'G') == ref)]]
  
  # conditional prob
  
  
  PA = PA(pr)
  PT = PT(pr)
  PC = PC(pr)
  PG = PG(pr)
  
  p.cond = c(dmultinom(c(Watson, Crick), prob = c(PA, PT)/2), # A
             dmultinom(c(Watson, Crick), prob = c(PT, PA)/2), # T
             dmultinom(c(Watson, Crick), prob = c(PC, PG)/2), # C
             dmultinom(c(Watson, Crick), prob = c(PG, PC)/2), # G
             dmultinom(c(Watson, Crick), prob = c(PA+PC, PT+PG)/4), # AC
             dmultinom(c(Watson, Crick), prob = c(PA+PG, PT+PC)/4), # AG
             dmultinom(c(Watson, Crick), prob = c(PA+PT, PT+PA)/4), # AT
             dmultinom(c(Watson, Crick), prob = c(PC+PG, PG+PC)/4), # CG
             dmultinom(c(Watson, Crick), prob = c(PC+PT, PG+PA)/4), # CT
             dmultinom(c(Watson, Crick), prob = c(PG+PT, PC+PA)/4)  # GT
  )
  
  # posterior prob
  p.post.unnorm = p.cond*theta
  p.post = p.post.unnorm/sum(p.post.unnorm)
  
  # prob not mutation (same with ref)
  # regarded as p.value
  
  p.value = p.post[1:4][which.max(c('A', 'T', 'C', 'G') == ref)]
  return(c(p.post, p.value, sum(Watson), sum(Crick), sum(p.post[1:4])))
}

# test
prob.post = matrix(0, nrow = nrow(d), ncol = length(STATUS) + 4)
# status.pred = rep('N', nrow(d))


for (i in 1:nrow(d)) {
  pp = postp(d[[i,2]], d[i,4] == 'CG', Watson[i,], Crick[i,])
  # pp = postp(ref = d[i,2], cg = d[i,4] == 'CG', c(0,0,2,2), c(0,0,2,2))
  
  # names(pp) = STATUS
  # barplot(pp, ylim = c(0,1))
  # 
  prob.post[i,] = pp
  # status.pred[i] = STATUS[which.max(pp)]
  
}


# max.prob = rowMaxs(prob.post[,1:10])
# sum(max.prob<0.95)


## allele frequencies

allele.weights = t(matrix(
  c(1, 0, 0, 0, 0.5, 0.5, 0.5, 0  , 0  , 0  ,
    0, 1, 0, 0, 0,   0  , 0.5, 0  , 0.5, 0.5,
    0, 0, 1, 0, 0.5, 0  , 0  , 0.5, 0.5, 0  ,
    0, 0, 0, 1, 0  , 0.5, 0  , 0.5, 0  , 0.5
  ),
  nrow = 4, byrow = T
))

allele.freq = prob.post[,1:10] %*% allele.weights
