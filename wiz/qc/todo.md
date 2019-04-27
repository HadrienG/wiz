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

result coding region
    min max 0m6,936s result ok  31,396s
    region  0m7.002s result +1  31,973s

Possibility of sequence with a unique contig, dendrotetra bug if a only one sequence in contig

[X] change implementation of hmmdb,
[X] write class contigs,
[X] cleaning code,
[X] comment code,
[X] make test,
[X] add data in report,
