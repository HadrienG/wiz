# todo for folder qc
* [X] Moving graph function into a report.py file
* [X] Create branch tetra to develop the tetra module

## file bins.py :
None


## file qc.py :
* [X] fixing genome importing file tools.py
* [X] adding window size parameter
* [X] adding auto filter
* [X] adding ignorefilter parameter
* [X] adding auto import folder


## file gc.py :
### extract_value
* [X] Fix the new bug to extract_value
* [X] Fix the filtering function


## file tetra.py :
* [X] Give relative tetra distribution
* [X] calculate metric with euclidian or Manathan distant between some contigs (Chebyshev distance added to compare with the others)
* [ ] Enable the choice of distance calculating method between two sequences (by parameter)
* [X] implementing a multiprocessing method
* [ ] adding args to enable choice of cpu nbr used in this module


## file tools.py :
* [X] bugs when we import file 
* [X] try to import file by file
* [X] \-> try to import group of files # * OK >> Bug here 
* [X] \-> try to import folder of files # * OK >> Bug when we import a folder with a unique file
* [X] \-> check if a gen is not already loaded # * OK
* [ ] \-> Check if file is fasta file # ! important way to not crash the code


## file report.py :
* [ ] Make repport in HTML with jinja2
### Distplot graph :
* [X] Adding information about gc distplot graph
* [X] Fixing the range of the yaxis 
* [X] Fixing ploting bug
### scatterplot graph :
* [X] Changing scatterplot graph in plot graph
* [ ] make a colorgroup with each couple sequence bounds
* [ ] mask name of bound in legend
### tetraplot graph :
* [X] Make a heatmap with proportion in each tetranuc
* [X] Make a heatmap with relative proportion for each tetranucl for each contig # ! Amount of data too large to display numeric values
* [ ] annotate the heatmap and specify the distance calculation method used
* [ ] sort sequence with hclust method



## to devtest task :
* [X] Split the E.coli genome in 100 part to test graphs (same process with S.cerevisiae)
* [X] Make a fake seq fasta to test data representation (A shuffle with the splited sequence from E.coli and S.cerevisiae)


## the reminders for the man with the head in the moon
joblib, principal composant analysis, MDS 



metric notes :
time pipenv run wiz qc -g /stage/Shuffle_200 -w 500  
    scipy funct    
    euclidian distance 26.890s   3st
    manhattan distance 26.525s   1st
    chebyshev distance 27.306s   4st
    minkowski distance 26.675s   2st

    without other func                      HomeMade Function
    euclidian distance 18.385s   3st            17.621s
    manhattan distance 18.257s   1st            17.294s
    chebyshev distance 18.972s   4st            
    minkowski distance 18.815s   2st

metrics notes about the func tetra.tetranuc_count and tetra.tetranuc_count_list
time pipenv run wiz qc -w 500 -g ../Shuffle_200 --debug
time to run the 2 func in a unique run : 1m32.994
                0                      : 0m01.187
    tetra.tetranuc_count(bin.seq)      : 0m06.911  likeable approach
    tetra.tetranuc_count_list(bin.seq) : 1m25.789  bad approach

metrics notes about the func tetra.distance_calculation and tetra.distance_MPcalculation
time pipenv run wiz qc -w 500 -g ../Shuffle_200 --debug
                            real        user
distance_MPcalculation      0m11s260    0m34s923
distance_calculation        0m24s586    0m24s595


result coding region    true result = 81,278007584
82.7624095139607
86.31937263012755
81.13796966563254
87.60556704584626
83.9861254739745
88.89391589107204
87.55386073767666
88.38331609789728
86.4701826956222
88.19588073078248
78.56127197518097
88.21742502585316
85.93803860737677
84.68846949327819
76.49301964839711
88.21527059634609
91.44045156842469
90.30075835918649
85.72905894519131
81.99112375043089
87.25654946570148
85.00086177180283
83.97319889693209
87.45044812133746
86.1168562564633

83.93872802481903
88.05153395380904
83.12866253016202
88.24974146845915
88.32083764219234
89.90003447087211
89.41528783178214
90.74672526714926
87.41813167873147
88.66554636332299
86.62745604963806
89.35927266459841
87.53877973112719
86.69208893485005
77.41511547742158
89.78369527749052
91.96182350913477
92.47026887280249
87.56247845570493
84.91683902102723
87.97182006204757
85.56101344364012
87.92011375387797
88.6590830748018
89.35927266459841