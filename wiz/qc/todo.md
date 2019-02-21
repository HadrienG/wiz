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
* [ ] Enable the choix of distance calculating method between two sequences (by parameter)


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
* [ ] make a group with each sequence bounds
### tetraplot graph :
* [X] Make a heatmap with proportion in each tetranuc
* [ ] Make a chart with relative proportion for each tetranucl for each contig # ! Amount of data too large to display numeric values
* [ ] annotate the heatmap and specify the distance calculation method used



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
