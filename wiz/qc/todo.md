# todo for folder qc
	*[X] Moving graph function into a report.py file
	*[ ] Create branch tetra to develop the tetra module

## file bins.py :
None


## file qc.py :
	[X] fixing genome importing file tools.py
	[X] adding window size parameter
	[X] adding auto filter
	[X] adding ignorefilter parameter
	[X] adding auto import folder
	[ ] adding some parameters


## file gc.py :
	### extract_value ###
	[X] Fix the new bug to extract_value
	[X] Fix the filtering function


## file tetra.py :
	[ ] Give relative tetra distribution
	[ ] Search metric with euclidian or Manathan distant between some contigs

## file tools.py :
	[X] bugs when we import file 
	* [X] try to import file by file
	* [X] \-> try to import group of files # * OK >> Bug here 
	* [X] \-> try to import folder of files # * OK >> Bug when we import a folder with a unique file
	* [X] \-> check if a gen is not already loaded # * OK
	* [ ] \-> Check if file is fasta file # ! important way to not crash the code

## file report.py :
[ ] Make repport in HTML with jinja2
	### Distplot graph :
	* [X] Adding information about gc distplot graph
	* [X] Fixing the range of the yaxis 
	* [X] Fixing ploting bug
	### scatterplot graph :
	* [X] Changing scatterplot graph in plot graph
	### tetraplot graph :
	* [ ] Make a barchart with proportion in each tetranuc
	* [ ] Make a chart with relative proportion for each tetranucl for each contig

