# Network Analysis with Cytoscape via py4cytoscape
These two scripts demonstrate how to automate network analysis with Cytoscape via the Python client py4cytoscape.

The EBSeq script demonstrates how to perform network analysis with a gene list after DE inference, querying the STRING database using given DE genes alone.

The rMATS script demonstrates how to focus on a particular disease, download genes related to a particular disease from the DISEASE database, run network analysis then highlight the relevant genes from an external gene list (which in this case are rMATS results from a splicing study).


## Notes
Cytoscape is a GUI and requires X11 forwarding if it is to be run on an HPC. The scripts also require the Cytoscape apps stringApp and MCODE which needed to be installed separately.

Cytoscape automatically logs progress in the user home directory, which could quickly overload it. Regular clean-up or re-direction of the logs are needed.

From experience, Cytoscape struggles with gene lists longer than 1000 genes (and the networks will be visually messy anyway), so large gene lists need to be filtered first.

![Up-regulated](/img/up-regulated.jpg)

![Down-regulated](/img/down-regulated.jpg)

![Splicing](/img/splicing.jpg)
