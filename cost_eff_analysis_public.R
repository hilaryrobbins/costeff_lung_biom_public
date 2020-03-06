
##############
# Title: Cost-effectiveness of biomarker-informed eligibility for lung cancer screening
# Code author: Hilary A. Robbins (International Agency for Research on Cancer, Lyon, France)
# Contact: RobbinsH@iarc.fr
# Start date: 16 January 2019
# Last edit date: 6 March 2020
##############

# READ ME! 

# How to use this code:
# This code corresponds to the paper "Assessment of Biomarker Testing for Lung Cancer Screening Eligibility"
  # JAMA Netw Open. 2020;3(3):e200409. doi:10.1001/jamanetworkopen.2020.0409
  # https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2762044
# By running this code in R, you can change the assumptions made in the primary analysis
  # and examine the impact they have on the results and conclusions.
# Throughout, the abbreviation fp/FP is false-positive, and tp/TP is true-positive. 
# The "risk.based" scenario is the reference scenario, which assumes screening eligibility based on a risk model accounting for smoking and demographic information.
# How to change key assumptions? Look for the objects "number.tested.with.biomarker", "biomarker.costs", and "cost.workup" (cost.workup can be calculated using different assumptions about false-positives.)

rm(list=ls(all=TRUE))  #Clears console
packages <- c("tidyverse","gmodels","scales")
lapply(packages, require, c = T)
options(scipen=999999)

#### The first section below stores and/or creates some key parameters underlying the analysis. ####

# First, save a few fixed parameters related to the US National Lung Screening Trial (NLST team NEJM 2011 and Black NEJM 2014)
people.screened <- 26642  # Black p1795
avg.screens.per.person <- (26309+24715+24102)/26722  # NLST NEJM 2011 Table 2
screens <- people.screened*avg.screens.per.person # estimated for Black population
perc.screens.as.fp <- 17497/(26309+24715+24102)  # NLST NEJM 2011 Table 3 (numerator) and Table 2 (denominator)
fp.screens <- perc.screens.as.fp*screens # estimated for Black population
all.cases <- 1076 # Black p1796
tp.screens <- 649 # NLST NEJM 2011 Table 3 - we assume same for Black
# The following are direct costs. Screening costs are per person, workup costs are calculated per positive screen, and treatment costs are calculated per case.
cost.screening.per.person <- 1130 # Black Table 2 
cost.workup.per.positive <- (835*people.screened)/(fp.screens+tp.screens) # Black Table 2
cost.treatment.per.case <- (1106*people.screened)/all.cases # Black Table 2. Note - we are omitting $3/person for radiation induced cancers
# The following are life-years in the presence and absence of screening
ly.case.screen <- 8.4792  # Black Table 1 lifetime horizon
ly.case.no.screen <- 6.8479  # Black Table 1 lifetime horizon
ly.noncase.screen <- 15.0097  # Black Table 1 lifetime horizon
ly.noncase.no.screen <- 15.0103  # Black Table 1 lifetime horizon

# Save a few parameters related to the USA population of ever-smokers and the number screening eligible (from Katki JAMA 2014)
num.to.screen <- 9018693 # Table 4 - using 3rd column
num.us.smokers <- 43413257 # Table 4
avg.lcrat.risk.all.smokers <- .015 # Table 4
avg.lcrat.risk.uspstf <- .039 # Table 4
avg.lcrat.risk.risk.based <- .045 # Table 4

# Calculated parameters
avg.lcrat.risk.nlst <- .0279  # Mean LCRAT risk in the CT arm of NLST (calculated by H. Robbins in analysis_nlst_split_v24_falsepos_only_v3.R)
cases.per.person <- all.cases/people.screened
sensitivity <- tp.screens/all.cases
perc.pos.as.tp <- tp.screens/(fp.screens+tp.screens)
est.cases.risk.based <- (avg.lcrat.risk.risk.based/avg.lcrat.risk.nlst)*cases.per.person*num.to.screen
all.cases.smokers <- (avg.lcrat.risk.all.smokers/avg.lcrat.risk.nlst)*cases.per.person*num.us.smokers
perc.all.cases.risk.based <- est.cases.risk.based/all.cases.smokers

# How many ever-smokers undergo biomarker testing to ultimately select 9.02 million?
number.tested.with.biomarker <- c(num.us.smokers/2)  # assume we test half of ever-smokers - this is a key assumption that can be varied

#### The next section performs the bulk of the calculations using the parameters generated above. ####

# Create overall dataframe structure (scenarios)
smooth.perc.all.cases = c(perc.all.cases.risk.based+.0001, perc.all.cases.risk.based+.0005, round(seq(perc.all.cases.risk.based+.001, 0.801, .001), 3))
cedat <- data.frame(
  scenario = c("all.smokers","uspstf","risk.based",rep("add", times=length(smooth.perc.all.cases))),
  percent.all.cases = c(NA, NA, NA, smooth.perc.all.cases))
# Add parameters to these scenarios
cedat <- mutate(cedat,
                avg.lcrat.risk = c(avg.lcrat.risk.all.smokers, avg.lcrat.risk.uspstf, avg.lcrat.risk.risk.based, rep(NA, nrow(cedat)-3)),
                people.screened = c(num.us.smokers, rep(num.to.screen, nrow(cedat)-1)),
                total.screens = people.screened*avg.screens.per.person,
                num.cases = ifelse(scenario %in% c("all.smokers","uspstf","risk.based"), # In these 3 scenarios, use average risk to scale up/down the NLST cases.per.person
                                   (avg.lcrat.risk/avg.lcrat.risk.nlst)*cases.per.person*people.screened,
                                    percent.all.cases*all.cases.smokers),  # Else, calculate based on the hypothetical percentage of all cases
                percent.all.cases = ifelse(scenario %in% c("all.smokers","uspstf","risk.based"), # For the first 3 scenarios, calculate the percentage of all cases captured
                                   (num.cases/all.cases.smokers), 
                                   percent.all.cases),                     # Else, leave the hypothetical quantity
                true.positives = sensitivity*num.cases,  # Calculate TPs by assuming a constant percentage of cases would be detected
                false.positives.lower = perc.screens.as.fp*total.screens,  # Calculate lower bound for FPs by assuming constant rate of false-pos/all-screens (not true - Kovalchik NEJM)
                false.positives.upper = (true.positives/perc.pos.as.tp) - true.positives, # Calculate upper bound for FPs by assuming constant rate of true-pos/all-pos (not true - Kovalchik NEJM)
                false.positives.mean = (false.positives.lower+false.positives.upper)/2)  # For the calculation below, we'll use the mean of the upper and lower bounds.
# Make a subset dataset that only includes the hypothetical scenarios
cedat.hyp <- dplyr::filter(cedat, substr(scenario, 1, 3)=="add")

# cedat.sub forms the core data frame of hypothetical scenarios. Stack it on top of itself to create scenarios for the biomarker costs and number of ever-smokers tested.
biomarker.costs <- c(5,50,100,150,200,250,300)  # list of possible biomarker costs - this can be changed to examine the impact of different per-test biomarker costs.
cedat.exp1 <- bind_rows(replicate(length(biomarker.costs), cedat.hyp, simplify=F))
cedat.exp1 <- mutate(cedat.exp1, cost.biomarker.per.test = rep(biomarker.costs, each=nrow(cedat.hyp)))
cedat.exp2 <- bind_rows(replicate(length(number.tested.with.biomarker), cedat.exp1, simplify=F))
cedat.exp2 <- mutate(cedat.exp2, number.tested = rep(number.tested.with.biomarker, each=nrow(cedat.exp1)))

# Add the 3 fixed scenarios to this dataframe
cedat.fix <- dplyr::filter(cedat, (scenario %in% c("all.smokers","uspstf","risk.based"))) %>% mutate(cost.biomarker.per.test=0, number.tested=0)
cedat.exp2 <- bind_rows(cedat.fix, cedat.exp2)

# Calculate the costs and life-years for each scenario
cedat.exp2 <- mutate(cedat.exp2,
                     cost.screening = cost.screening.per.person*people.screened,
                     cost.workup = cost.workup.per.positive*(true.positives+false.positives.mean),  # To assume different numbers (extremes/bounds) of false-positives, can substitue .mean with .lower or .upper (see above)
                     cost.treatment = cost.treatment.per.case*num.cases,
                     cost.biomarker = cost.biomarker.per.test*number.tested,
                     cost.total = cost.screening + cost.workup + cost.treatment + cost.biomarker,
                     ly.screen = ly.case.screen*num.cases + ly.noncase.screen*(people.screened - num.cases),
                     ly.no.screen = ly.case.no.screen*num.cases + ly.noncase.no.screen*(people.screened - num.cases),
                     ly.gained = ly.screen - ly.no.screen) 

# Save cost and ly.gained from risk-based scenario and use to calculate ICERs
ref.ly.gained <- (cedat.exp2 %>% dplyr::filter(scenario=="risk.based") %>% select(ly.gained))[1,1]
ref.cost.total <- (cedat.exp2 %>% dplyr::filter(scenario=="risk.based") %>% select(cost.total))[1,1]

# Calculate ICERs
cedat.exp2 <- mutate(cedat.exp2,
                     incremental.cost = ifelse(scenario=="uspstf", NA, cost.total - ref.cost.total),
                     incremental.ly.gained = ifelse(scenario=="uspstf", NA, ly.gained - ref.ly.gained),
                     icer = ifelse(scenario=="uspstf", NA, incremental.cost/incremental.ly.gained))  # calculate ICER as the incremental cost per incremental LYG


#### You can now examine your results in any way you like - they are stored in the dataframe cedat.exp2. ####
# For example:
cedat.exp2 %>% filter(percent.all.cases==0.7)
cedat.exp2 %>% filter(percent.all.cases==0.65)
cedat.exp2 %>% filter(percent.all.cases==0.8)

#### The following code plots the ICERs in the same manner as Figure 1 in the paper. ####

top.plot <- 150000  # change this to change the cap on the y-axis for the plot.
cedat.exp2.plot <- cedat.exp2 %>% dplyr::filter(substr(scenario, 1, 3)=="add")  # include only the hypothetical biomarker scenarios in plot data
dpi <- 300  # set resolution
png(file="file_path_here", width=16.7*dpi, height=11.8*dpi, res=dpi)
ggplot(data=cedat.exp2.plot, aes(x=percent.all.cases, y=icer)) + geom_line(aes(group=cost.biomarker.per.test)) + theme_bw() +
  theme(axis.title = element_text(size=30), axis.text = element_text(size=22), panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  coord_cartesian(ylim=c(10000, top.plot)) +
  scale_x_continuous(labels = scales::percent_format(accuracy=1), breaks=seq(0.50,0.80,0.05)) + 
  scale_y_continuous(labels = dollar, breaks=seq(25000,top.plot,25000)) + 
  geom_label(data=dplyr::filter(cedat.exp2.plot, percent.all.cases==.66), aes(label = paste("$",cost.biomarker.per.test,sep="")), size=8) +
  xlab("\nPercentage of ever-smoking lung cancer cases captured into a\nfixed-size screening population using different eligibility strategies") +
  ylab("Incremental cost-effectiveness ratio (ICER)\n") + 
  geom_text(aes(x=0.715, y=135500, label="< Per-person biomarker cost for\n biomarker-informed strategies"), size=8) +
  geom_hline(yintercept=50000, linetype="dashed", size=0.8) +
  geom_segment(linetype="solid", size=0.3, y=0, yend=170000, x=cedat.exp2[cedat.exp2$scenario=="risk.based","percent.all.cases"], xend=cedat.exp2[cedat.exp2$scenario=="risk.based","percent.all.cases"]) +
  annotate(geom="text", hjust=0, x=0.545, y=125000, size=6.5, label="< USPSTF eligibility\n54% of cases") +
  geom_segment(linetype="solid", size=0.3, y=0, yend=170000, x=cedat.exp2[cedat.exp2$scenario=="uspstf","percent.all.cases"], xend=cedat.exp2[cedat.exp2$scenario=="uspstf","percent.all.cases"]) +
  annotate(geom="text", hjust=1, x=0.62, y=100000, size=6.5, label="Smoking-model eligibility >\n(reference strategy)\n62% of cases") +
  NULL
dev.off()

