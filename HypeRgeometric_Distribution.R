###########################################################################################
#Name: Jay S. Gandhi                                                                      #
#Advisors: Dr. Bradford Windle, PhD & Ms. Tara Nulton MS                                  #
#Purpose: R script to perform the Hypergeometric Distribution to calculate the            #
#         rates of co-mutation of patient samples with TP53 (non-silent) gene mutations   # 
#         in Lung Squamous Cell Carcinoma (TCGA-LUSC) and Lung Adenocarcinoma (TCGA-LUAD).# 
###########################################################################################

#Begin by importing the tables in to the environment

#Create an object for p53 table
message("Choose p53 file")
p53_list <- read.csv(file.choose(), header = TRUE)

#Create an object for Sample table
message("Choose the Sample file")
sample_list <- read.csv(file.choose(), header = TRUE)

#loop to calculate the total count of p53 sample ids
p_total = 0
for (i in 1:length(p53_list)){
  p_col_count = sum(p53_list[,i]!="")
  p_total = p_col_count + p_total
}

#loop to calculate the total count of sample ids
s_total = 0
for (s in 1:length(sample_list)){
  s_col_count = sum(sample_list[,s]!="")
  s_total = s_col_count + s_total
}

#declare a new data frame to store the outputs
hyper.df = data.frame()

#total patient id counts (i.e. p53 sample ids + sample ids)
#total_count = 1054
#update 01/09/2018
#NOTE: 406 is the total P53 sample IDs after removing 9 repeative sample IDs and redudant sample IDs
total_count = 1027

#for loop to traverse each column of the table
for (gene in 1:length(sample_list)){

  #get the name of the gene (header) from the sample
  gene_name = (names(sample_list[gene]))
  
  #get the length of the gene sample list
  sample_count = sum(sample_list[,gene]!="")
  
  #for loop to traverse each column of the TP53 clusters
  for (p in 1:length(p53_list)){

    #get the name of the p53 header from the p53 table
    p53_cluster_name = (names(p53_list[p]))
    
    #get the number of patient ids in the p53 cluster
    p53_cluster_count = sum(p53_list[,p]!="")
    
    #calculate the expected value 
    #NOTE: 433 is the total P53 sample IDs
    
    #update 01/09/2019
    #NOTE: 406 is the total P53 sample IDs after removing 9 repeative sample IDs and redudant sample IDs
    expected_intersect = (((p53_cluster_count)/(406))*(sample_count))

    #intersect each sample with an individual cluster
    observed_intersect = (intersect(p53_list[,p],sample_list[,gene]))
    #obtains the count of the sample IDs without the whitespace
    observed_intersect = length(observed_intersect[observed_intersect != ""])
    
    #conditional statement to calculate either a left tail or a right tail hypergeometric distribution
    #if the observed (intersect value) is less than the expected
    if (observed_intersect < expected_intersect){
      
      #left tail phyper calculation
      p_val = phyper(observed_intersect, sample_count, total_count - sample_count, p53_cluster_count, lower.tail = TRUE)    
    
    } else {
      
      #right tail phyper calculation 
      p_val = phyper(observed_intersect - 1, sample_count, total_count - sample_count, p53_cluster_count, lower.tail = FALSE)
    }
    
    #bind newly calculated values to columns
    new.hyper.df = cbind(gene_name, p53_cluster_name, total_count, observed_intersect, expected_intersect, sample_count, p53_cluster_count, p_val)
    
    #bind the one row of data to previous data frame
    hyper.df = rbind(hyper.df, new.hyper.df)
    
  }
}
#change column names
names(hyper.df) = c("Gene Name", "TP53 Cluster Name", "Total Samples", "Observed Intersect", "Expected Intersect", "Sample Count", "TP53 Cluster Count", "P-Value")

#export the dataset to an output file
write.csv(hyper.df, "TCGA_LUSC_LUAD_Hypergeometric_Output.csv")
